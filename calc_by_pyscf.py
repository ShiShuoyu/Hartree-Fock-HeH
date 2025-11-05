from pyscf import gto

# 直接计算，指定为闭壳层
mol = gto.Mole()
mol.unit = 'AU'
mol.atom = 'H 0 0 0; He 0 0 1.4632'
mol.basis = 'STO-3G'
mol.spin = 0 # 闭壳层
mol.charge = 1
mol.build()

overlap_matrix = mol.intor('int1e_ovlp')
print(overlap_matrix)
