a_Newton = (inner(grad(v), dot(as_tensor(dj_ddgT), grad(T_)))+
            inner(grad(v), as_vector(dj_ddT))*T_)*dxm
res = -inner(grad(v), j)*dxm
