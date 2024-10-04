#utilisation de la méthode de Jacobi pour la résolution d'un système linéaire
#' methode de jacobi
#'
#'La methode de jacobi, due au mathematicien allemand Karl Jacobi, est une methode iterative
#'de resolution d un systeme matriciel de la forme Ax=B. Pour cela, on utilise une suite x^k
#'qui converge vers un point fixe x, solution du systeme d equations lineaires.
#'On cherche a construire l algorithme pour x^0 donne, la suite x^(k+1)=F(x^k) avec k appartenant N.
#'
#' @author DJEUMENI Alexia
#' @param matrice une matrice carree
#' @param xo approximation initiale
#' @param B second membre du systeme
#' @param nomb_ite nombres d iterations a effectuer
#'
#' @examples jacobi(matrix(c(9,1,1,4,6,-2,1,0,-6),3,3),
#' matrix(c(0,0,0),3,1),matrix(c(-17,4,14),3,1),2)



jacobi <- function(matrice,xo,B,nomb_ite){ #création d'une fonction qui va
  #nous permettre de réaliser aisément nos calculs
  #fonction qui permet de créer la matrice inférieure de la matrice initiale
  #avec les éléments de la diagonale nulles
  infe <- function(matr){
    for(i in 1:nrow(matr))
      for (j in 1:ncol(matr)) {
        if(i<j){
          matr[i,j] <- 0
        }
      }
    diag(matr) <- 0
    return(matr)
  }
  #fonction qui permet de créer la matrice supérieure de la matrice initiale
  #avec les éléments de la diagonale nulles
  sup <- function(matr){
    for(i in 1:nrow(matr))
      for (j in 1:ncol(matr)) {
        if(i>j){
          matr[i,j] <- 0
        }
      }
    diag(matr) <- 0
    return(matr)
  }
  #ici on montre  la convergence de la matice
  converge <- TRUE
  if(det(matrice)!=0 & nrow(matrice) == nrow(xo) & ncol(xo)==1 & nrow(matrice)==nrow(B) & ncol(B)==1 ){ #conditions pour qu'on puisse vérifier la convergence de la matrice
    for (i in 1:nrow(matrice)) {
      condition <- 0
      lin <- matrice[i,]
      k <- abs(matrice[i,i])
      for (j in lin) {
        condition <- condition + abs(j)
      }
      condition <- condition - k
      if(k<condition){
        converge <- FALSE
      }
    }
    if(converge==FALSE){
      #s'il n'y a pas convergence,  le message ci-dessous s'affichera
      return("pas de convergence")
    }
    else{
      #dans le cas contraire, c'est le message ci-dessous qui va s'afficher
      D <- diag(diag(matrice))
      E <- infe(matrice)
      F <- sup(matrice)
      #la boucle for va permettre de calculer les itérés de Jacobi
      for (i in 1:nomb_ite) {
        xo <- solve(D)%*%((-1)*(E+F)%*%xo + B)
      }

    }

  }
  else
    return("erreur de synthaxe")
  return(xo)
}




#jacobi(matrix(c(9,1,1,4,6,-2,1,0,-6),3,3),matrix(c(0,0,0),3,1),matrix(c(-17,4,14),3,1),2)
#jacobi(matrix(c(3,1,1,-2),2,2),matrix(c(0,0),2,1),matrix(c(-1,4),2,1),3)
#jacobi(matrix(c(9,1,1,4,-2,6,1,0,-6),3,3),matrix(c(0,0,0),3,1),matrix(c(-17,4,14),3,1),4)
