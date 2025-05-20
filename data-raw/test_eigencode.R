library(Matrix)

L1=sparseMatrix(c(1,2),c(3,3),dims = c(3,3))
L1=-(L1+t(L1))
diag(L1)=c(1,1,2)
L1
L2=sparseMatrix(c(1,3),c(2,4),dims = c(4,4))
L2=-(L2+t(L2))
diag(L2)=rep(1,4)
L2


cutset=rbind(c(0,3),c(1,6))
perm=0:6
L=gtclust:::buildLaplacianR(L1,L2,cutset,perm)

L11=gtclust:::remove1row1col(L)
Matrix::determinant(L[-1,-1],log=TRUE)$modulus

ch=Matrix::Cholesky(L11,LDL = TRUE)
sum(log(diag(as(ch,"triangularMatrix"))))

lu=qr(L[,-1])


# 2  0 -1 -1  0  0  0
# 0  2 -1  0  0  0 -1
# -1 -1  2  0  0  0  0
# -1  0  0  2 -1  0  0
# 0  0  0 -1  1  0  0
# 0  0  0  0  0  1 -1
# 0 -1  0  0  0 -1  2
# Log nb tree
# 12.7526



# 1  0 -1  0  0  0  0  0  0
# 0  3 -1  0 -1 -1  0  0  0
# -1 -1  2  0  0  0  0  0  0
# 0  0  0  2  0  0 -1 -1  0
# 0 -1  0  0  2 -1  0  0  0
# 0 -1  0  0 -1  3 -1  0  0
# 0  0  0 -1  0 -1  3 -1  0
# 0  0  0 -1  0  0 -1  3 -1
# 0  0  0  0  0  0  0 -1  1
# Log nb tree
# -2

ij=rbind(c(1,3),
c(2,3),
c(2,5),
c(2,6),
c(4,7),
c(4,8),
c(5,6),
c(6,7),
c(7,8),
c(8,9))
L=sparseMatrix(ij[,1],ij[,2],dims = c(9,9))
L=(L+t(L))
dd=rowSums(L)
L=-L
diag(L)=dd
L
L11=gtclust:::remove1row1col(L)
Matrix::determinant(L11,log=TRUE)$modulus
ch=Matrix::Cholesky(L11,LDL = TRUE)
sum(log(diag(as(ch,"triangularMatrix"))))

# 1  0 -1  0  0  0  0  0  0  0  0  0
# 0  2 -1  0  0  0  0  0  0  0  0 -1
# -1 -1  2  0  0  0  0  0  0  0  0  0
# 0  0  0  1  0 -1  0  0  0  0  0  0
# 0  0  0  0  3 -1  0 -1 -1  0  0  0
# 0  0  0 -1 -1  2  0  0  0  0  0  0
# 0  0  0  0  0  0  2  0  0 -1 -1  0
# 0  0  0  0 -1  0  0  2 -1  0  0  0
# 0  0  0  0 -1  0  0 -1  3 -1  0  0
# 0  0  0  0  0  0 -1  0 -1  3 -1  0
# 0  0  0  0  0  0 -1  0  0 -1  3 -1
# 0 -1  0  0  0  0  0  0  0  0 -1  2
# Log nb tree
# 29.2212

ij=rbind(c(1,3),
         c(2,3),
         c(2,12),
         c(4,6),
         c(5,6),
         c(5,8),
         c(5,9),
         c(7,10),
         c(7,11),
         c(8,9),
         c(9,10),
         c(10,11),
         c(11,12))
L=sparseMatrix(ij[,1],ij[,2],dims = c(12,12))
L=(L+t(L))
dd=rowSums(L)
L=-L
diag(L)=dd
L
Matrix::determinant(L[-1,-1],log=TRUE)$modulus



L1 = L[1:6,1:6]
L2= L[7:12,7:12]


chol1=chol(L1[-1,-1])
chol2=chol(L2[-1,-1])
