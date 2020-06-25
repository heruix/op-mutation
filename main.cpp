#include <stdio.h>
#include <stdlib.h>
#include <cstdint>
#include <future>
#include <sstream>

#include "flint/flint.h"
#include "flint/fmpq_mat.h"
#include <set>
#include <random>
#include <utility>

#include <fstream>
#include <iostream>

#include "arith.h"

#define sqr(x) ((x)*(x))

constexpr uint64_t test_cases[] = { 
	0, 1
};
constexpr size_t tc_size = sizeof(test_cases) / sizeof(uint64_t);

std::random_device device{};

OpNode* gen_node(int depth)
{
	if (depth > 0)
	{
		const auto new_binop = new BinopNode();
		new_binop->type = static_cast<BinopNode::BinopType>(std::uniform_int_distribution<int>(0, static_cast<int>(BinopNode::BinopType::BINOP_MAX) - 1)(device));
		new_binop->right = std::shared_ptr<OpNode>(gen_node(depth - 1));
		new_binop->left = std::shared_ptr<OpNode>(gen_node(depth - 1));
		return new_binop;
	}
	else
	{
		switch (std::uniform_int_distribution<int>(0, 2)(device))
		{
			case 0:
			{
				const auto new_unop = new UnopNode();
				new_unop->type = static_cast<UnopNode::UnopType>(std::uniform_int_distribution<int>(0, static_cast<int>(UnopNode::UnopType::UNOP_MAX) - 1)(device));
				new_unop->right = std::static_pointer_cast<OpNode>(std::make_shared<VarNode>(std::uniform_int_distribution<int>(0, 1)(device)));
				return new_unop;
			}
			case 1:
			{
				return new VarNode(std::uniform_int_distribution<int>(0, 1)(device));
			}
		} 
	}
}

void gen_nodes(int depth, unsigned count, std::vector<std::shared_ptr<OpNode>>& nodes)
{
	for (unsigned i = 0; i < count; i++)
		nodes.push_back(std::shared_ptr<OpNode>(gen_node(depth)));
}

void rank_factorization(fmpz_mat_t A, fmpz_mat_t Ar, int rank, fmpz_mat_t C, fmpz_mat_t F) {
	//C = pivot columns of A
	int l = 0;
	auto b = std::make_unique<bool[]>(Ar->c);
	for (int i = 0; i < Ar->c; i++) {
		if (l < Ar->r && !fmpz_is_zero(fmpz_mat_entry(Ar, l, i))) {
			b[i] = true;
			l++;
		}
	}

	fmpz_mat_init(C, Ar->r, l);
	for (int i = 0, j = 0; i < Ar->c; i++) 
	{
		if (!b[i])
			continue;
		for (int k = 0; k < Ar->r; k++)
			fmpz_set(fmpz_mat_entry(C, k, j), fmpz_mat_entry(A, k, i));
		j++;
	}

	//F = nonzero rows of Ar
	fmpz_mat_init(F, rank, Ar->c);
	for (int k = 0; k < rank; k++)
		for (int j = 0; j < Ar->c; j++)
			fmpz_set(fmpz_mat_entry(F, k, j), fmpz_mat_entry(Ar, k, j));
}


static int64_t f64(uint64_t x, uint64_t y)
{
	return y * (2 - y * x);
}

static uint64_t findInverse64(uint64_t x)
{
	uint64_t y = (3 * x) ^ 2;
	y = f64(x, y);
	y = f64(x, y);
	y = f64(x, y);
	y = f64(x, y);
	return y;
}

bool gen_variable(std::function<uint64_t(uint64_t, uint64_t)> test, std::stringstream& out, int recurse = 0) {
	fmpz_mat_t A, B, AB{};

	constexpr auto m_size = sqr(tc_size);
	constexpr auto fun_size = 50;

	std::vector<std::shared_ptr<OpNode>> nodes{};
	gen_nodes(1, fun_size, nodes);

	fmpz_mat_init(A, m_size, fun_size);
	fmpz_mat_init(B, m_size, 1);

	fmpz_mat_init(AB, m_size, fun_size + 1);
	
	for (unsigned i = 0; i < tc_size; i++)
		for (unsigned j = 0; j < tc_size; j++)
			for (unsigned k = 0; k < fun_size; k++) {
				const auto a = nodes[k]->compute(test_cases[i], test_cases[j]);
				fmpz_set_ui(fmpz_mat_entry(A, i * tc_size + j, k), a);
				fmpz_set_ui(fmpz_mat_entry(AB, i * tc_size + j, k), a);
			}

	for (unsigned i = 0; i < tc_size; i++)
		for (unsigned j = 0; j < tc_size; j++) {
			const auto a = test(test_cases[i], test_cases[j]);
			fmpz_set_ui(fmpz_mat_entry(B, i * tc_size + j, 0), a);
			fmpz_set_ui(fmpz_mat_entry(AB, i * tc_size + j, fun_size) , a);
		}

	fmpz_mat_t Ar;
	fmpz_mat_init_set(Ar, A);

	fmpz_t mmod;
	fmpz_init_set_ui(mmod, 0x100000000);
	fmpz_mul_ui(mmod, mmod, 0x100000000);

	fmpz_t den;
	fmpz_init(den);
	auto rank = fmpz_mat_rref(Ar, den, A);

	if (rank != fmpz_mat_rank(AB)) {
		fmpz_mat_clear(A);
		fmpz_mat_clear(B);
		fmpz_mat_clear(AB);
		return false;
	}
	fmpz_mat_t API;
	fmpz_mat_init(API, A->c, A->r);

	//moore-penrose pseudoinverse
	//sources: wikipedia
	
	fmpz_mat_t AT, AM;
	fmpz_mat_init(AT, A->c, A->r);
	fmpz_mat_transpose(AT, A);
	fmpz_mat_init(AM, A->r, AT->c);
		
	if (rank == A->r) {
		fmpz_t den;
		fmpz_mat_mul(AM, A, AT);
		fmpz_mat_inv(AM, den, AM);

		fmpz_mat_mul(API, AT, AM);
	} 
	else if (rank == A->c) {
		fmpz_t den;
		fmpz_mat_mul(AM, AT, A);
		fmpz_mat_inv(AM, den, AM);

		fmpz_mat_mul(API, AM, AT);	
	} 
	else {
		fmpz_mat_t C, F;
		rank_factorization(A, Ar, rank, C, F);

		//scuffed moore-penrose pseudoinverse thing
		//thanks wikipedia
		
		//rank decompose A -> A = CF
		//pinv(A) = (inv(C^T * C) * C^T) * (F^T * inv(F * C^T))
		//lol
		
		//CT = C^T
		//FT = F^T
		fmpz_mat_t CT, FT;

		fmpz_mat_init(CT, C->c, C->r);
		fmpz_mat_transpose(CT, C);

		fmpz_mat_init(FT, F->c, F->r);
		fmpz_mat_transpose(FT, F);

		//CM = C * CT, CI = inv(CM)
		//FM = F * FT, FI = inv(FM)
		fmpz_mat_t CM, CI, FM, FI;
		fmpz_t den;

		fmpz_mat_init(CM, CT->r, C->c);
		fmpz_mat_init(CI, CT->r, C->c);

		fmpz_mat_mul(CM, CT, C);
		fmpz_mat_inv(CI, den, CM);
		
		fmpz_mat_init(FM, F->r, FT->c);
		fmpz_mat_init(FI, F->r, FT->c);

		fmpz_mat_mul(FM, F, FT);
		fmpz_mat_inv(FI, den, FM);

		//CM2, FM2 = CI * CT, FT * FI
		fmpz_mat_t CM2, FM2;
		fmpz_mat_init(CM2, CI->r, CT->c);
		fmpz_mat_init(FM2, FT->r, FI->c);
		
		fmpz_mat_mul(CM2, CI, CT);
		fmpz_mat_mul(FM2, FT, FI);

		fmpz_mat_scalar_mod_fmpz(CM2, CM2, mmod);
		fmpz_mat_scalar_mod_fmpz(FM2, FM2, mmod);
			
		//API = FM2 * CM2
		fmpz_mat_mul(API, FM2, CM2);

		fmpz_mat_clear(CM2);
		fmpz_mat_clear(FM2);

		fmpz_mat_clear(CT);
		fmpz_mat_clear(CM);
		fmpz_mat_clear(CI);

		fmpz_mat_clear(FT);
		fmpz_mat_clear(FM);
		fmpz_mat_clear(FI);
	}

	fmpz_mat_t T;
	fmpz_mat_init(T, A->r, API->c);
	fmpz_mat_mul(T, A, API);
	fmpz_mat_scalar_mod_fmpz(T, T, mmod);

	mp_limb_t nmod_offs = 1;
	auto r = &T->rows[0][0];
	uint64_t rv = fmpz_get_ui(r);
	if (!fmpz_is_one(r))
	{
		if (rv % 2 == 1)
		{
			nmod_offs = findInverse64(rv);
		}
		else
		{
			return false;
		}
	}

	fmpz_mat_t X;
	fmpz_mat_init(X, API->r, B->c);
	fmpz_mat_mul(X, API, B);
	fmpz_mat_scalar_mod_fmpz(X, X, mmod);
	
	out << "(";
	bool first = true;
	for (int j = 0; j < fun_size; j++) {
		const uint64_t variable = fmpz_get_ui(fmpz_mat_entry(X, j, 0));
		if (variable == 0)
			continue;

		if (first)
			first = false;
		else
			out << " + ";

		out << "(" << variable << "ULL * " << nodes[j]->to_string() << ")";
	}
	out << ") * " << nmod_offs << "ULL";

	fmpz_mat_clear(A);
	fmpz_mat_clear(B);
	fmpz_mat_clear(AB);

	fmpz_mat_clear(X);

	fmpz_mat_clear(AT);
	fmpz_mat_clear(AM);

	return true;
}

int main(int argc, char** argv)
{
	printf("OPERATION HEADASSIFIER\n");
	const auto test = [](uint64_t a, uint64_t b) -> uint64_t {
		return a + b;
	};

	//printf("%d\n", g2(15, 15));

	while (true) {
		std::stringstream out;
		while (!gen_variable(test, out, 0)) { }
		printf("%s\n\n", out.str().c_str());
	}
}	