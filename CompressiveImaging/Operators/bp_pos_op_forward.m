function OUT = bp_pos_op_forward(x, A, W, N)

OUT = [A(W(x(1:N))); W(x(1:N)) - x(N+1:2*N)];
