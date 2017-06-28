function OUT = bp_pos_op_backward(x, At, Wt, M)

OUT = [Wt(At(x(1:M))) + Wt(x(M+1:end)); - x(M+1:end)];