#pragma once

#include <memory>
#include <string>

struct OpNode
{
	std::shared_ptr<OpNode> left;
	std::shared_ptr<OpNode> right;

	virtual std::string to_string() { return "[NO]"; };
	virtual uint64_t compute(uint64_t a, uint64_t b) { return -1; };
	virtual ~OpNode() = default;
};

struct BinopNode : public OpNode
{
	enum BinopType
	{
		AND,
		OR,
		XOR,
		ADD,
		BINOP_MAX,
		MUL,
	};

	BinopType type;
	std::string to_string() override
	{
		std::string op;
		switch (type) {
		case ADD:
			op = "+";
			break;
		case MUL:
			op = "*";
			break;
		case AND:
			op = "&";
			break;
		case OR:
			op = "|";
			break;
		case XOR:
			op = "^";
			break;
		//case MOD:
		//	return "(" + left->to_string() + " % (" + right->to_string() + " | 1))";
		default:
			op = "UNK";
			break;
		}
		return "(" + left->to_string() + " "  + op + " " + right->to_string() + ")";
	}

	uint64_t compute(uint64_t a, uint64_t b) override
	{
		uint64_t l = left->compute(a, b);
		uint64_t r = right->compute(a, b);
		switch (type) {
		case ADD:
			return l + r;
		case MUL:
			return l * r;
		case AND:
			return l & r;
		case OR:
			return l | r;
		case XOR:
			return l ^ r;
		//case MOD:
		//	return l % (r | 1);
		default:
			return -1;
		}
	}
};

struct UnopNode : public OpNode
{
	enum UnopType
	{
		NOT,
		UNOP_MAX,
		MINUS,
		NONE,
	};

	UnopType type;
	std::string to_string() override
	{
		std::string op;
		switch (type) {
		case MINUS:
			op = "-";
			break;
		case NOT:
			op = "~";
			break;
		default:
			op = "";
			break;
		}
		return op + right->to_string();
	}

	uint64_t compute(uint64_t a, uint64_t b) override
	{
		uint64_t r = right->compute(a, b);
		switch (type) {
		case MINUS:
			return static_cast<uint64_t>(-static_cast<int32_t>(r));
		case NOT:
			return ~r;
		default:
			return r;
		}
	}
};

struct VarNode : public OpNode
{
	bool isA;
	std::string to_string() override
	{
		return isA ? "a" : "b";
	}

	uint64_t compute(uint64_t a, uint64_t b) override
	{
		return isA ? a : b;
	}

	VarNode(bool a) : isA(a) {}
};

struct ConstantNode : public OpNode
{
	uint64_t constant = 1;// __rdtsc() % static_cast<uint64_t>(-1);
	std::string to_string() override
	{
		return std::to_string(constant) + "UL";
	}

	uint64_t compute(uint64_t a, uint64_t b) override
	{
		return constant;
	}
};