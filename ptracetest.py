
def ptrace(psi,sel):
    import numpy as np
    drho = psi.dims[0]
    M = 1
    N = 1
    prod_drho_sel = 1
    prod_drho_rest = 1
    cnt = 0
    rest = np.zeros(len(drho)-len(sel))
    for i in range(len(drho)):
        N = N*drho[i]
        insel = False
        for j in range(len(sel)):
            if (i==sel[j]): 
                M=M*drho[i]
                prod_drho_sel = prod_drho_sel*drho[i]
                insel = True
        if (not insel):
            prod_drho_rest = prod_drho_rest*drho[i]
            rest[cnt] = i
            cnt = cnt+1

    #for i in range(len(sel)):
    #    psi_permuted[i] = psi[sel[i]]
    #for i in range(len(sel),len(psi)):
    #    psi_permuted[i] = psi[rest[i-len(sel)]]


    print prod_drho_sel
    print prod_drho_rest
    print psi.data.data
    A = np.reshape(psi.data.toarray(),[prod_drho_rest,prod_drho_sel]).T
    return dot(conjg(A.T),A)

