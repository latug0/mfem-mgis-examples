a_Newton = -(v*dH1_ddT*T_-
             dt*theta*(inner(grad(v), dot(as_tensor(dj1_ddgT), grad(T_)))+
                       inner(grad(v), as_vector(dj1_ddT))*T_))*dxm
res = (v*(H1-H0)-dt*(theta*inner(grad(v), j1)+(1-theta)*inner(grad(v), j0)))*dxm