function y_i = rk4(dgl, t, y, u, h, para)

    k1 = dgl(t,     y,          u, para);
    k2 = dgl(t+h/2, y+0.5*k1*h, u, para);
    k3 = dgl(t+h/2, y+0.5*k2*h, u, para);
    k4 = dgl(t+h  , y+    k3*h, u, para);
    y_i = y + h/6*(k1+2*k2+2*k3+k4);

end