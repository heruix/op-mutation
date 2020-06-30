# op-mutation
experiments in mixed-boolean arithmetic obfuscation

# what is mixed-boolean arithmetic?
a term some people smarter than me coined for arithmetic expressions that contain both traditional boolean algebra (and, or, not) and arithmetic (add, sub, mul) operators

# how does this work?
quite a bit of math. the truth tables of random operations are placed as column vectors in a matrix. the moore-penrose pseudoinverse is then used to compute the solution to the linear equation Ax = B, where A is the random operation matrix and B is the truth table of the target operation. 
if the pseudoinverse is accurate, multiplying the original matrix by the psuedoinverse will lead to an identity matrix. however, if the matrix does not have a proper inverse, it will lead to diagonal matrices with coefficients that aren't one. 
in this case, we can check for if the coefficient is odd. if it isn't, there's nothing we can do and we have to bail. otherwise, we can multiply the entire statement by the modular inverse of the coefficient in the matrix, allowing us to continue.
after that, the inverted matrix is multiplied with B to find the solution vector x. in this case, x represents a vector of coefficients, one for each random operation in the matrix. if all is well, the resulting equation looks something like this:
```(a1 * c1 + a2 * c2 + a3 * c3 ... ) * modular_inverse```
where a1, a2, a3, ... correspond to the values in the solution vector, c1, c2, c3, ... correspond to the random operations, and modular_inverse is the previously explained modular inverse, or 1 if the inverse is proper.

# boring, show me something funny
how about a [4.5 megabyte addition](examples/add)

# additional information
due to recent improvements in how https://github.com/vtil-project/VTIL-Core handles arithmetic expressions in multiplication, this mutation generator is almost completely broken by it :) 
