# Function to calculate matrix and constrains for quadratic programing
# Non-overlapping strata
PSWquadra.matrix = function(w, p.hat, X, Z, epsilon, a, b, TwoSideBound = FALSE) {
  require(quadprog)
  n = length(Z)
  n0 = sum(Z==0)
  n1 = sum(Z==1)
  epsilon = (n0*n1*epsilon)/(n0+n1)
  if (is.null(w)) {
    Dmat = matrix(0, n, n)
    diag(Dmat) = 1
  } else if(!is.null(w)) {
    w = matrix(w, ncol = 1)
    Dmat = w%*%w
  }
  w.hat = 1/(Z*p.hat + (1-Z)*p.hat)
  #w.hat = weight.calculate.ps(Z, p.hat, standardize = TRUE)$weight
  dvec = -w.hat
  I0 = which(Z == 0)
  I1 = which(Z == 1)
  A.eq = rep(1,n)
  A.eq[I0] = 1
  A.eq[I1] = -1
  XZ = X*(Z-(1-Z))
  # Box Bound
  a = bound.jit(a)
  b = bound.jit(b)
  abound = Z/b + (1-Z)/(1-a)
  bbound = Z/a + (1-Z)/(1-b)
  if (TwoSideBound == FALSE) {
    Amat = cbind(A.eq, XZ, -XZ, diag(n), -diag(n))
    bvec = c(0,rep(-epsilon, ncol(XZ)*2), abound, -bbound)
  } else if (TwoSideBound == TRUE) {
    require(Matrix)
    Amat = cbind(A.eq, XZ,diag(n))
    Amat = t(Amat)
    #Dmat = t(Dmat)
    Dmat = Matrix(Dmat, sparse = TRUE)
    Amat = Matrix(Amat, sparse = TRUE)
    l = c(0,rep(-epsilon,ncol(XZ)),abound)
    u = c(0,rep(epsilon,ncol(XZ)),bbound)
    bvec = list(l = l, u = u)
  }
  out = list(Dmat = Dmat, dvec = dvec, Amat = Amat,
             bvec = bvec)
  return(out)
}

bound.jit = function(a) {
  if (a == 0) {
    a = a + 1e-06
  } else if (a == 1) {
    a = a - 1e-06
  } else {
    a = a
  }
  return(a)
}
# Function to estimate the weight with non-overlap strata
PSWquadra.nonoverlap.est = function(p.hat, X, Z, elist = list(input = input, e.min = e.min, bylength = bylength, a = a),
                                    ej, sol_name) {
  data = cbind(p.hat, Z)
  estEps = MaxEpsilonFind(X = X, Z = Z, p.hat = p.hat, ej = ej, a = elist$a)
  strataIndex = list()
  X.fit = list()
  epsilon = list()
  for(i in 1:nrow(ej)) {
    I = which(p.hat >= ej[i,1] & p.hat < ej[i,2])
    strataIndex[[i]] = I
    if (length(I) == 1) {
      #X0 = standardize(matrix(X[I,],nrow = 1))$X
      X0 = matrix(X[I,],nrow = 1)
    } else {
      X0 = standardize(X[I,])$X
    }
    X.fit[[i]] = cbind(data[I,], X0)
    if (is.null(elist$input)) {
      e.max = log(estEps[i,3],2)
      e.temp = c(seq(elist$e.min, e.max, elist$bylength), e.max)
      e.temp = unique(e.temp)
      epsilon[[i]] = 2^e.temp
    } else {
      epsilon[[i]] = elist$input
    }
  }
  X.out = cbind(data,X)
  w = list()
  e.out = c()
  for (i in 1:length(X.fit)) {
    p = length(strataIndex[[i]])
    if (p > 0 && p < 5) {
      w0 = rep(1, nrow(X.fit[[i]]))
      e0 = NA
    } else {
      for (j in 1:length(epsilon[[i]])) {
        print(paste0("i = ",i,"; epsilon = ", epsilon[[i]][j]))
        if (sol_name == "quadprog") {
          quadmatrix = PSWquadra.matrix(w = NULL, p.hat = X.fit[[i]][,1],
                                        X = X.fit[[i]][,3:ncol(X.fit[[i]])],
                                        Z = X.fit[[i]][,2], epsilon = epsilon[[i]][j],
                                        a = ej[i,1], b = ej[i,2],
                                        TwoSideBound = FALSE)
          quad.fit = NULL
          tryCatch({quad.fit = solve.QP(Dmat = quadmatrix$Dmat, dvec = quadmatrix$dvec,
                                        Amat = quadmatrix$Amat, bvec = quadmatrix$bvec, meq = 1)},
                   error = function(e) {
                     quad.fit = NULL
                   })
          if (is.null(quad.fit)) {
            w0 = rep(NA, nrow(X.fit[[i]]))
            e0 = epsilon[[i]][j]
          } else if (!is.null(quad.fit)) {
            w0 = quad.fit$solution
            e0 = epsilon[[i]][j]
            break
          }
        } else if (sol_name == "osqp") {
          quadmatrix = PSWquadra.matrix(w = NULL, p.hat = X.fit[[i]][,1],
                                        X = X.fit[[i]][,3:ncol(X.fit[[i]])],
                                        Z = X.fit[[i]][,2], epsilon = epsilon[[i]][j],
                                        a = ej[i,1], b = ej[i,2],
                                        TwoSideBound = TRUE)
          quad.fit = NULL
          tryCatch({quad.fit = osqp(P = quadmatrix$Dmat, q = quadmatrix$dvec,
                                    A = quadmatrix$Amat, l = quadmatrix$bvec$l,
                                    u = quadmatrix$bvec$u,
                                    pars = osqpSettings(max_iter = 10000L))$Solve()},
                   error = function(e) {
                     quad.fit = NULL
                   })
          if (is.null(quad.fit)) {
            w0 = rep(NA, nrow(X.fit[[i]]))
            e0 = epsilon[[i]][j]
          } else if (!is.null(quad.fit) && quad.fit$info$status_val != 1) {
            w0 = rep(NA, nrow(X.fit[[i]]))
            e0 = epsilon[[i]][j]
          } else if (!is.null(quad.fit) && quad.fit$info$status_val == 1) {
            w0 = quad.fit$x
            e0 = epsilon[[i]][j]
            break
          }
        }
      }
    }
    w[[i]] = w0
    e.out[i] = e0
  }
  e.out = cbind(ej, e.out)
  colnames(e.out) = c("a","b","epsilon")
  out = list(w = w, e.out = e.out, epsilon = epsilon, strataIndex = strataIndex)
  return(out)
}

# Main function to call in non-overlap strata estimation
PSWquadra.fit = function(p.hat, X, Z, epsilon = list(input = NULL, e.min = -10, bylength = 0.4, a = 3),
                         ej, sol = list(sol_name = "osqp", min.sub = 0.1)) {
  require(osqp)
  n = length(Z)
  sample.size = sample.size.calculate(Z = Z, p.hat = p.hat, ej = ej)
  n.sub.min = ceiling(n*sol$min.sub)
  In = unique(which(sample.size[,1] <= 5 | sample.size[,3] <= 5 | sample.size[,4] <= 5))
  if (length(In) == 0) {
    ej.fit = ej
    error = c(0,0)
    no.sub = 0
  } else {
    no.sub = sum(sample.size[In,1])
    if (no.sub <= n.sub.min) {
      ej.fit = ej[-In,]
      error = c(0,no.sub)
    } else if (no.sub > n.sub.min) {
      ej.fit = ej[-In,]
      error = c(1,no.sub)
    }
  }
  fit = PSWquadra.nonoverlap.est(p.hat = p.hat, X = X, Z = Z, elist = epsilon,
                                 ej = ej.fit, sol_name = sol$sol_name)
  strataIndex = fit$strataIndex
  success0 = c()
  w.hat = weight.calculate.ps(Z = Z, ps = p.hat, standardize = FALSE)$weight
  w.out = w.hat
  for(i in 1:nrow(ej.fit)) {
    w.out[strataIndex[[i]]] = fit$w[[i]]
    if (sum(is.na(fit$w[[i]])) != 0) {
      success0[i] = 0
    } else if (sum(is.na(fit$w[[i]])) == 0) {
      success0[i] = 1
    }
  }
  success = rep(NA, nrow(ej))
  if(length(In) == 0) {
    success = success0
  } else {
    success[-In] = success0
  }
  sample.size = cbind(sample.size, "success" = success)
  #I.hybrid = which(sample.size[,3] <= 10 | sample.size[,4] <= 10 & sample.size[,5] == 0)
  I.hybrid = which(sample.size[,5] == 0)
  no.sub.hyb = sum(sample.size[I.hybrid,1])
  if (length(I.hybrid) == 0) {
    w.out.hybrid = w.out
    no.sub.all = 0 + no.sub
  } else if (length(I.hybrid) != 0 && no.sub.hyb > (n.sub.min-no.sub)) {
    w.out.hybrid = w.out
    no.sub.all = 0 + no.sub
    error[1] = 1
  } else if (length(I.hybrid) == 1 && no.sub.hyb <= (n.sub.min-no.sub)) {
    I2 = which(p.hat >= ej[I.hybrid,1] & p.hat < ej[I.hybrid,2])
    w.out.hybrid = w.out
    w.out.hybrid[I2] = w.hat[I2]
    no.sub.all = length(I2) + no.sub
  } else if (length(I.hybrid) != 1 && no.sub.hyb <= (n.sub.min-no.sub)) {
    I2 = which(p.hat >= ej[I.hybrid[1],1] & p.hat < ej[I.hybrid[1],2])
    for(i in 2:length(I.hybrid)) {
      I2 = c(I2, which(p.hat >= ej[I.hybrid[i],1] & p.hat < ej[I.hybrid[i],2]))
    }
    I2 = unique(I2)
    w.out.hybrid = w.out
    w.out.hybrid[I2] = w.hat[I2]
    no.sub.all = length(I2) + no.sub
  }
  e.out = cbind(ej, rep(NA,nrow(ej)))
  if(length(In) == 0) {
    e.out[,3] = fit$e.out[,3]
  } else {
    e.out[-In,3] = fit$e.out[,3]
  }
  error = c(error, no.sub.all)
  names(error) = c("error","no.sub","no.sub.hybrid")
  out = list(weight = w.out, weight.hybrid = w.out.hybrid, sample.size = sample.size,
             epsilon = e.out, p.hat = p.hat, epsilon.list = fit$epsilon, error = error)
  return(out)
}
# Overlap strata
PSWquadra.matrix.OverlapStrata = function(w = NULL, p.hat, X, Z, epsilon, ej, X_global,
                                          TwoSideBound = FALSE) {
  require(osqp)
  require(Matrix)
  require(quadprog)
  n = length(Z)
  #n_effect = ceiling(max(EffectiveSampleSize(Z = Z, p.hat = p.hat, ej = ej)))
  n_effect = EffectiveSampleSize(Z = Z, p.hat = p.hat, ej = ej)
  if (is.null(w)) {
    Dmat = matrix(0, n, n)
    diag(Dmat) = 1
  } else if(!is.null(w)) {
    w = matrix(w, ncol = 1)
    Dmat = w%*%w
  }
  w.hat = 1/(Z*p.hat + (1-Z)*p.hat)
  dvec = -w.hat
  # Generate new covariate matrix with local information
  k = nrow(ej)
  Kij = matrix(rep(0,n*k), nrow = n, ncol = k)
  for(i in 1:k) {
    Kij[which(p.hat >= ej[i,1] & p.hat < ej[i,2]),i] = 1
  }
  #X.fit = standardize(X*Kij[,1])$X
  X.fit = XKij_standardize(X = X, Kij = Kij, a = 1)
  for(i in 2:k) {
    #X0 = standardize(X*Kij[,i])$X
    X0 = XKij_standardize(X = X, Kij = Kij, a = i)
    X.fit = cbind(X.fit, X0)
  }
  intercept.fit = Kij
  I0 = which(Z == 0)
  I1 = which(Z == 1)
  A.eq = rep(1,n)
  A.eq[I0] = 1
  A.eq[I1] = -1
  A.eq = A.eq*intercept.fit
  if (is.null(X_global)) {
    XZ = X.fit*(Z-(1-Z))
  } else {
    XZ_global = X_global*(Z-(1-Z))
    XZ_local = X.fit*(Z-(1-Z))
    XZ = cbind(XZ_local, XZ_global)
  }
  # Generate box constraints
  abound.matrix = matrix(nrow = n, ncol = k)
  bbound.matrix = matrix(nrow = n, ncol = k)
  for(i in 1:n) {
    for(j in 1:k) {
      a = bound.jit(ej[j,1])
      b = bound.jit(ej[j,2])
      if(intercept.fit[i,j]==0) {
        abound.matrix[i,j] = NA
        bbound.matrix[i,j] = NA
      } else if(intercept.fit[i,j]==1) {
        abound.matrix[i,j] = Z[i]/b + (1-Z[i])/(1-a)
        bbound.matrix[i,j] = Z[i]/a + (1-Z[i])/(1-b)
      }
    }
  }
  abound = c()
  bbound = c()
  for(i in 1:n) {
    abound[i] = max(abound.matrix[i,],na.rm = T)
    bbound[i] = min(bbound.matrix[i,],na.rm = T)
  }
  if (is.null(X_global)) {
    n_effect_vec = rep(n_effect, each = ncol(X))
  } else {
    n_e_all = (as.numeric(table(Z)[1])*as.numeric(table(Z)[2]))/length(Z)
    n_effect_vec = c(rep(n_effect, each = ncol(X)),rep(n_e_all, ncol(X_global)))
  }
  if (length(epsilon) == 1) {
    epsilon.vec0 = rep(-epsilon, ncol(XZ))
    epsilon.vec = rep(-epsilon, ncol(XZ)*2)
  } else if (length(epsilon) == k) {
    r = ncol(X)
    epsilon.vec0 = rep(epsilon[1],r)
    for(i in 2:k) {
      epsilon.vec0 = c(epsilon.vec0, rep(epsilon[i],r))
    }
    epsilon.vec0 = -epsilon.vec0
    #epsilon.vec = rep(epsilon.vec0,2)
    if (is.null(X_global)) {
      epsilon.vec0 = epsilon.vec0
    } else {
      e.m = median(epsilon, na.rm = T)
      epsilon.vec0 = c(epsilon.vec0, rep(-e.m,ncol(X_global)))
    }
    epsilon.vec = rep(epsilon.vec0,2)
  }
  epsilon.vec = epsilon.vec*c(n_effect_vec,n_effect_vec)
  epsilon.vec0 = epsilon.vec0*n_effect_vec
  #epsilon.vec = epsilon.vec*n_effect
  #epsilon.vec0 = epsilon.vec0*n_effect
  if (TwoSideBound == FALSE) {
    Amat = cbind(A.eq, XZ, -XZ, diag(n), -diag(n))
    bvec = c(rep(0,k), epsilon.vec, abound, -bbound)
  } else if (TwoSideBound == TRUE) {
    Dmat = Matrix(Dmat, sparse = TRUE)
    Amat = cbind(A.eq, XZ,diag(n))
    Amat = t(Amat)
    Amat = Matrix(Amat, sparse = TRUE)
    l = c(rep(0,k),epsilon.vec0,abound)
    u = c(rep(0,k),-epsilon.vec0,bbound)
    bvec = list(l = l, u = u)
  }
  out = list(Dmat = Dmat, dvec = dvec, Amat = Amat,
             bvec = bvec, k = k)
  return(out)
}

OverlapStrata.estimate = function(w= NULL, p.hat, X, Z, epsilon, ej, X_global,
                                  sol_name) {
  if (sol_name == "quadprog") {
    quadmatrix = PSWquadra.matrix.OverlapStrata(w = NULL, p.hat, X, Z, epsilon, ej, X_global, TwoSideBound = FALSE)
    quad.fit = NULL
    tryCatch({quad.fit = solve.QP(Dmat = quadmatrix$Dmat, dvec = quadmatrix$dvec,
                                  Amat = quadmatrix$Amat, bvec = quadmatrix$bvec, meq = quadmatrix$k)},
             error = function(e) {
               quad.fit = NULL
             })
    if (is.null(quad.fit)) {
      w = rep(NA, nrow(X))
      e0 = epsilon
      #e.out = cbind(ej, e0, rep(i,k))
    } else if (!is.null(quad.fit)) {
      w = quad.fit$solution
      e0 = epsilon
    }
  } else if (sol_name == "osqp") {
    quadmatrix = PSWquadra.matrix.OverlapStrata(w = NULL, p.hat, X, Z, epsilon, ej, X_global, TwoSideBound = TRUE)
    quad.fit = NULL
    tryCatch({quad.fit = osqp(P = quadmatrix$Dmat, q = quadmatrix$dvec,
                              A = quadmatrix$Amat, l = quadmatrix$bvec$l,
                              u = quadmatrix$bvec$u,
                              pars = osqpSettings(max_iter = 10000L))$Solve()},
             error = function(e) {
               quad.fit = NULL
             })
    if (is.null(quad.fit)) {
      w = rep(NA, nrow(X))
      e0 = epsilon
    } else if (!is.null(quad.fit) && quad.fit$info$status_val != 1) {
      w = rep(NA, nrow(X))
      e0 = epsilon
    } else if (!is.null(quad.fit) && quad.fit$info$status_val == 1) {
      w = quad.fit$x
      e0 = epsilon
    }
  }
  out = list(w = w, e0 = e0)
}
# Main function to call in overlap strata estimation
PSWquadra.fit.OverlapStrata = function(p.hat, X, Z, epsilon = list(input = NULL, e.min = -10, bylength = 0.4, a = 3),
                                       ej, X_global, sol = list(sol_name = "osqp", epsilon.div = c(3,3,3), min.sub = 0.1)) {
  X = standardize(X)$X
  X.orig = X
  Z.orig = Z
  p.hat.orig = p.hat
  n = length(Z)
  sample.size = sample.size.calculate(Z = Z, p.hat = p.hat, ej = ej)
  In = unique(which(sample.size[,1] <= 5 | sample.size[,3] <= 5 | sample.size[,4] <= 5))
  if (length(In) == 0) {
    ej.fit = ej
    error = c(0,0)
    X = X.orig
    Z = Z.orig
    p.hat = p.hat.orig
    epsilon.div = sol$epsilon.div
  } else {
    n.sub.min = ceiling(n*sol$min.sub)
    ej.fit = ej[-In,]
    e.div = sol$epsilon.div
    e.div0 = rep(seq(1,length(e.div)), e.div)
    e.div0 = e.div0[-In]
    epsilon.div = as.numeric(table(e.div0))
    if (length(epsilon.div) == length(e.div)) {
      epsilon.div = epsilon.div
    } else if (length(epsilon.div) < length(e.div)) {
      r = sum(epsilon.div)/length(e.div)
      if (r <= 1) {
        epsilon.div = rep(1,sum(epsilon.div))
      } else {
        r1 = round(r)
        epsilon.div0 = rep(r1, (length(e.div) - 1))
        r2 = sum(epsilon.div) - sum(epsilon.div0)
        epsilon.div = c(epsilon.div0, r2)
        rm(epsilon.div0, r1, r2)
      }
    }
    I.include.exclude = sample.in.strata(p.hat = p.hat, ej = ej, In = In)
    no.sub = I.include.exclude$no.exclude
    if (no.sub != 0 && no.sub <= n.sub.min) {
      X = X.orig[I.include.exclude$sample.include, ]
      Z = Z.orig[I.include.exclude$sample.include]
      p.hat = p.hat.orig[I.include.exclude$sample.include]
      error = c(0,no.sub)
    } else if (no.sub > n.sub.min) {
      X = X.orig[I.include.exclude$sample.include, ]
      Z = Z.orig[I.include.exclude$sample.include]
      p.hat = p.hat.orig[I.include.exclude$sample.include]
      error = c(1,no.sub)
    } else if (no.sub == 0) {
      X = X.orig
      Z = Z.orig
      p.hat = p.hat.orig
      error = c(0,no.sub)
    }
  }
  names(error) = c("error","no.sub")
  #data = cbind(p.hat, Z, X)
  sigma.na = EpsilonTildeEstimate(X = X, Z = Z, p.hat = p.hat, weight = NULL, ej = ej.fit)$epsilon_tilde_avg
  I.s.na = which(is.na(colMeans(sigma.na)))
  if (length(I.s.na) == 0) {
    X_global = X_global
    X = X
  } else {
    X_global = cbind(X[,I.s.na], X_global)
    X = X[,-I.s.na]
  }
  elist = epsilon
  estEps = MaxEpsilonFind(X = X, Z = Z, p.hat = p.hat, ej = ej.fit, a = elist$a)
  if (!is.null(dim(elist$input))) {
    epsilon.matrix = elist$input
    for (i in 1:nrow(epsilon.matrix)) {
      print(paste0("i = ",i,"; epsilon = ", epsilon.matrix[i,]))
      fit = OverlapStrata.estimate(w= NULL, p.hat, X, Z, epsilon = epsilon.matrix[i,],
                                   ej.fit, X_global, sol_name = sol$sol_name)
      w = fit$w
      e0 = c(fit$e0, rep(i,nrow(ej.fit)))
      #e.out = cbind(ej, e0, rep(i,nrow(ej)))
      if (sum(is.na(w)) == 0) {
        break
      }
    }
  } else if (!is.null(elist$input) && is.null(dim(elist$input))) {
    fit = PSW.epsilon.selection(p.hat = p.hat, X = X, Z = Z,
                                epsilon = elist$input, ej = ej.fit, X_global,
                                sol_name = sol$sol_name, e.div = epsilon.div)
    w = fit$w
    e0 = fit$epsilon
    #e.out = cbind(ej, e0)
  } else if (is.null(elist$input)) {
    e.max = log(max(estEps[,3], na.rm = T), 2)
    e.temp = c(seq(elist$e.min, e.max, elist$bylength), e.max)
    e.temp = unique(e.temp)
    epsilon = 2^e.temp
    fit = PSW.epsilon.selection(p.hat = p.hat, X = X, Z = Z,
                                epsilon = epsilon, ej = ej.fit, X_global = X_global,
                                sol_name = sol$sol_name, e.div = epsilon.div)
    w = fit$w
    e0 = fit$epsilon
    #e.out = cbind(ej, e0)
  }
  weight.hat = weight.calculate.ps(Z = Z, ps = p.hat.orig, standardize = FALSE)$weight
  if (length(In) == 0) {
    w.out = w
    e.out = cbind(ej, e0)
  } else {
    if (no.sub == 0) {
      w.out = w
      e.out = cbind(rep(NA, nrow(ej)), rep(NA, nrow(ej)))
      e.out[-In,] = e0
      e.out = cbind(ej, e.out)
    } else if (no.sub == n) {
      w.out = weight.hat
      e.out = cbind(ej, rep(NA, nrow(ej)), rep(NA, nrow(ej)))
    } else if (no.sub > 0 && no.sub < n) {
      w.out = weight.hat
      w.out[I.include.exclude$sample.include] = w
      e.out = cbind(rep(NA, nrow(ej)), rep(NA, nrow(ej)))
      e.out[-In,] = e0
      e.out = cbind(ej, e.out)
    }
  }
  colnames(e.out) = c("a","b","epsilon", "index")
  out = list(weight = w.out, sample.size = sample.size,
             epsilon = e.out, p.hat = p.hat.orig, epsilon.list = epsilon, error = error)
  return(out)
}

PSW.epsilon.selection = function(p.hat, X, Z, epsilon, ej, X_global, sol_name, e.div) {
  n = nrow(X)
  k = nrow(ej)
  for (i in 1:length(epsilon)) {
    print(paste0("i = ",i,"; epsilon = ", epsilon[i]))
    epsilon0 = rep(epsilon[i], k)
    fit = OverlapStrata.estimate(p.hat = p.hat, X = X, Z = Z, epsilon = epsilon0,
                                 ej = ej, X_global = X_global, sol_name = sol_name)
    w = fit$w
    e0 = fit$e0
    I = i
    if (sum(is.na(fit$w)) == 0) {
      break
    }
  }
  epsilon.out = cbind(e0, rep(I,k))
  e0_uni = unique(e0)
  if (e0_uni == epsilon[1]) {
    out = list(w = w, epsilon = epsilon.out)
  } else if (e0_uni == epsilon[length(epsilon)] && sum(is.na(w)) == n) {
    out = list(w = w, epsilon = epsilon.out)
  } else {
    e.new = epsilon[1:I]
    out = tree.selection(p.hat, X, Z, e.new, ej, X_global, w.init = w,
                         sol_name = sol_name, e.div = e.div)
  }
  return(out)
}

Index.generate = function(I.start, k) {
  I = apply(matrix(I.start, nrow = 1), 2, rep, each = k)
  I0 = diag(k)
  I1 = I - I0
  I1.na = ifelse(I1 == 0, 0, 1)
  I.na = which(rowMeans(I1.na) == 1)
  I1 = I1[I.na,]
  return(I1)
}

Index.generate.unique = function(I, k) {
  require(mgcv)
  I.out = Index.generate(I[1,],k)
  for(i in 2:nrow(I)) {
    I.out = rbind(I.out, Index.generate(I[i,],k))
  }
  out = uniquecombs(I.out)
  return(out)
}

tree.selection = function(p.hat, X, Z, e.new, ej, X_global, w.init, sol_name, e.div) {
  n = nrow(X)
  k0 = nrow(ej)
  k = length(e.div)
  p = length(e.new)
  min.sum = k - 1 + p
  I.start = Index.generate(rep(p, k),k)
  e.start = t(apply(I.start, 1, FUN = function(epsilon, div, x) rep(epsilon[x], times = div),
                    epsilon = e.new, div = e.div))
  epsilon.out0 = cbind(rep(e.new[p],k0), rep(p,k0))
  w = w.init
  while(TRUE) {
    w.new = list()
    success = c()
    for (i in 1:nrow(e.start)) {
      fit = OverlapStrata.estimate(p.hat = p.hat, X = X, Z = Z, epsilon = e.start[i,],
                                   ej = ej, X_global = X_global, sol_name = sol_name)
      w.new[[i]] = fit$w
      if (sum(is.na(fit$w)) == length(fit$w)) {
        success[i] = 0
      } else if (sum(is.na(fit$w)) == 0) {
        success[i] = 1
      }
    }
    if (sum(success) == 0) {
      w = w
      epsilon.out = epsilon.out0
      break
    } else {
      I.start0 = I.start[which(success == 1), ]
      if (sum(success) == 1) {
        w = w.new[[which(success == 1)]]
        epsilon.out0 = cbind(e.new[I.start0], I.start0)
        epsilon.out0 = apply(epsilon.out0, 2, FUN = function(x, div) rep(x, times = div),
                             div = e.div)
        I.start = Index.generate(I.start0, k)
        I.start.out = sum(I.start)
      } else {
        w = w.new[which(success == 1)]
        local.bal = c()
        for (i in 1:length(w)) {
          local.bal[i] = mean(abs(Fitted.strata.diff(X, Z, p.hat, w[[i]], ej)[-1,]), na.rm = T)
        }
        w = w[[which.min(local.bal)]]
        epsilon.out0 = cbind(e.new[I.start0[which.min(local.bal),]],
                             I.start0[which.min(local.bal),])
        epsilon.out0 = apply(epsilon.out0, 2, FUN = function(x, div) rep(x, times = div),
                             div = e.div)
        I.start = Index.generate.unique(I.start0, k)
        I.start.out = mean(rowSums(I.start))
      }
      if (I.start.out < min.sum) {
        break
      }
      #e.start = t(apply(I.start, 1, FUN = function(epsilon, x) epsilon[x], epsilon = e.new))
      e.start = t(apply(I.start, 1, FUN = function(epsilon, div, x) rep(epsilon[x], times = div),
                        epsilon = e.new, div = e.div))
    }
  }
  out = list(w = w, epsilon = epsilon.out0)
  return(out)
}

FindBound = function(ej) {
  ej = rbind(c(0,1),ej)
  out = matrix(nrow = nrow(ej), ncol = 4)
  for(i in 1:nrow(ej)) {
    a = bound.jit(ej[i,1])
    b = bound.jit(ej[i,2])
    Z = 0
    out[i,1] = Z/b + (1-Z)/(1-a)
    out[i,2] = Z/a + (1-Z)/(1-b)
    Z = 1
    out[i,3] = Z/b + (1-Z)/(1-a)
    out[i,4] = Z/a + (1-Z)/(1-b)
  }
  colnames(out) = c("T0.a","T0.b","T1.a","T1.b")
  return(out)
}

stack.matrix = function(fit, strata = "5strata") {
  M = matrix(nrow = length(fit), ncol = 8)
  if (strata == "5strata") {
    out = list(ej1 = M, ej2 = M, ej3 = M, ej4 = M, ej5 = M)
    k = 5
    ej = cbind(seq(0,0.8,0.2),seq(0.2,1,0.2))
  } else if (strata == "10strata") {
    out = list(ej1 = M, ej2 = M, ej3 = M, ej4 = M, ej5 = M,
               ej6 = M, ej7 = M, ej8 = M, ej9 = M, ej10 = M)
    k = 10
    ej = cbind(seq(0,0.9,0.1),seq(0.1,1,0.1))
  }
  for(i in 1:length(fit)) {
    for(j in 1:k) {
      if(is.na(fit[[i]])) {
        out[[j]][i,] = c(ej[j,], 1000, 0, rep(NA, 4))
      } else {
        out[[j]][i,] = fit[[i]][j,]
      }
    }
  }
  return(out)
}

test.sample.size = function(fit) {
  nf = length(which(fit[,4]==0))
  if (nf == 0) {
    sum.stat = c(NA, NA, mean(fit[,5],na.rm = T),NA,
                 sd(fit[,5],na.rm = T),NA,
                 NA, NA, NA,NA, NA,mean(fit[,6],na.rm = T),NA,
                 sd(fit[,6],na.rm = T),NA, 1,
                 mean(fit[,3], na.rm = T), sd(fit[,3], na.rm = T))
  } else if (nf == 1) {
    success = fit[which(fit[,4]==1),]
    fail = matrix(fit[which(fit[,4]==0),], nrow = 1)
    X2 = chisq.test(matrix(c(sum(fail[,7], na.rm = T), sum(fail[,8], na.rm = T),
                             sum(success[,7], na.rm = T), sum(success[,8], na.rm = T)),
                           nrow = 2, byrow = T))
    OR = X2$observed[1,1]*X2$observed[2,2]/(X2$observed[1,2]*X2$observed[2,1])
    sum.stat = c(NA, NA, mean(success[,5],na.rm = T),mean(fail[,5],na.rm = T),
                 sd(success[,5],na.rm = T), 0,
                 X2$statistic, X2$p.value, OR,
                 NA, NA, mean(success[,6],na.rm = T),mean(fail[,6],na.rm = T),
                 sd(success[,6],na.rm = T), 0, nrow(success)/nrow(fit),
                 mean(success[,3], na.rm = T), sd(success[,3], na.rm = T))
  } else {
    success = fit[which(fit[,4]==1),]
    fail = fit[which(fit[,4]==0),]
    ttest = t.test(success[,5], fail[,5])
    X2 = chisq.test(matrix(c(sum(fail[,7], na.rm = T), sum(fail[,8], na.rm = T),
                             sum(success[,7], na.rm = T), sum(success[,8], na.rm = T)),
                           nrow = 2, byrow = T))
    OR = X2$observed[1,1]*X2$observed[2,2]/(X2$observed[1,2]*X2$observed[2,1])
    prop.diff = t.test(success[,6], fail[,6])
    sum.stat = c(ttest$statistic, ttest$p.value, as.numeric(ttest$estimate),
                 sd(success[,5],na.rm = T),sd(fail[,5],na.rm = T),
                 X2$statistic, X2$p.value, OR,
                 prop.diff$statistic, prop.diff$p.value, as.numeric(prop.diff$estimate),
                 sd(success[,6],na.rm = T),sd(fail[,6],na.rm = T), nrow(success)/nrow(fit),
                 mean(success[,3], na.rm = T), sd(success[,3], na.rm = T))
  }
  names(sum.stat) = c("n.ttest.stat","n.ttest.pvalue","n.success.mean","n.fail.mean",
                      "n.success.sd","n.fail.sd","X2.stat","X2.pvalue","OR",
                      "prop.diff.stat","prop.diff.pvalue","case.success","case.fail",
                      "case.success.sd","case.fail.sd", "success.rate","success.epsilon.mean",
                      "success.epsilon.sd")
  return(sum.stat)
}

test.sample.size.allstrata = function(fit) {
  k = length(fit)
  N.col = paste0("ej",seq(1,k))
  out = matrix(nrow = 18, ncol = k)
  for(i in 1:k) {
    out[,i] = test.sample.size(fit[[i]])
  }
  colnames(out) = N.col
  rownames(out) = c("n.ttest.stat","n.ttest.pvalue","n.success.mean","n.fail.mean",
                    "n.success.sd","n.fail.sd","X2.stat","X2.pvalue","OR",
                    "prop.diff.stat","prop.diff.pvalue","case.success","case.fail",
                    "case.success.sd","case.fail.sd","success.rate","success.epsilon.mean",
                    "success.epsilon.sd")
  return(out)
}
#' Sample Size Calculation in Local Neighborhoods
#'
#' This function calculate the sample sizes in pre-specified local neighborhoods
#' (the number of subjects in each propensity score (PS) stratified sub-population).
#'
#' @param Z The binary treatment indicator. A vector with 2 unique numeric values in 0 = untreated and 1 = treated.
#' @param p.hat The propensity score used to determine the sub-population in each local neighborhood.
#' @param ej The matrix of the local neighborhoods, which contains two columns of positive values greater or equal to 0 and less or equal to 1.
#'           The rows of ej represent the neighborhoods. The first column is the start point of the local neighborhoods.
#'           The second column is the end point of the local neighborhoods.
#'
#' @return The matrix contains sample sizes in the local neighborhoods with each row corresponding to
#'         each local neighborhood (strata). The first gives the total number (n) of subjects in each strata.
#'         The third or fourth column give the sample sizes of untreated (n0) or treated (n1) subjects in each strata.
#'         The second column gives the ratio between n1 and n.
#' @examples
#' # Simulate data
#' KS = Kang_Schafer_Simulation(n = 1000, seeds = 5050)
#' # The treatment indicator and the true propensity score
#' Z = KS$Data[,2]
#' true.ps = KS$Data[,11]
#' # Local neighborhoods
#' ej = cbind(seq(0,0.7,0.1),seq(0.3,1,0.1))
#' # Calculate sample size in true PS-stratified sub-populations
#' local.sample.size = sample.size.calculate(Z = Z, p.hat = true.ps, ej)
#'
#' @export
sample.size.calculate = function(Z, p.hat, ej) {
  k = nrow(ej)
  out = matrix(nrow = k, ncol = 4)
  for(i in 1:k) {
    I = which(p.hat >= ej[i,1] & p.hat < ej[i,2])
    n_j = length(I)
    if (n_j == 0) {
      out[i,] = rep(0,4)
    } else if (n_j == 1) {
      Zj = Z[I]
      if (Zj == 0) {
        out[i,] = c(n_j, 0, 1, 0)
      } else if (Zj == 1) {
        out[i,] = c(n_j, 1, 0, 1)
      }
    } else {
      Zj = Z[I]
      if (sum(Zj) == length(Zj)) {
        out[i,] = c(n_j, 1, 0, length(Zj))
      } else if (sum(Zj) == 0) {
        out[i,] = c(n_j, 0, length(Zj), 0)
      } else {
        n0n1 = as.vector(table(Zj))
        out[i,] = c(n_j,n0n1[2]/length(Zj), n0n1)
      }
    }
  }
  colnames(out) = c("n","n1/n","n0","n1")
  return(out)
}

QuantileInMData = function(ps, k) {
  P = seq(0,1,1/k)
  P = P[-1]
  Q = matrix(nrow = ncol(ps), ncol = k)
  for(i in 1:ncol(ps)) {
    Q[i,] = quantile(ps[,i],probs = P)
  }
  out = colMeans(Q)
  #out[length(out)] = 1
  return(out)
}

BuildEpsilonMatrix = function(epsilon, ej) {
  k = nrow(ej)
  e_matrix = matrix(nrow = length(epsilon), ncol = k)
  for(i in 1:length(epsilon)) {
    e_matrix[i,] = rep(epsilon[i],k)
  }
  return(e_matrix)
}

#' Convert the Inverse Probability Weight (IPW) for ATE to Propensity Score
#'
#' This function converts the IPW weights from PSLB 2 to the corresponding propensity scores. When the subject i is
#' in the untreated group, the converted propensity score is 1-(1/weight_i). When the subject i is in the treated group,
#' converted propensity score is 1/weight_i.
#'
#' @param w The input IPW weights.
#' @param Z The binary treatment indicator. A vector with 2 unique numeric values in 0 = untreated and 1 = treated.
#'
#' @return The vector contanining the converted propensity score.
#'
#' @export
PS.ConvertFrom.Weight = function(w, Z) {
  probs.min = 1e-6
  w = as.vector(w)
  ps = rep(NA, length(w))
  ps[which(Z==1)] = 1/w[which(Z==1)]
  ps[which(Z==0)] = 1-(1/w[which(Z==0)])
  ps = pmin(1-probs.min,ps)
  ps = pmax(probs.min,ps)
  ps = as.vector(ps)
  return(ps)
}

SDInStrata = function(X, p.hat, ej) {
  sd = matrix(nrow = nrow(ej), ncol = ncol(X))
  for(i in 1:nrow(ej)) {
    I = which(p.hat >= ej[i,1] & p.hat < ej[i,2])
    if (length(I) == 1) {
      sd0 = rep(0,ncol(X))
    } else {
      X0 = X[I,]
      sd0 = apply(X0, 2, sd)
    }
    sd[i,] = sd0
  }
  return(sd)
}

XKij_standardize = function(X, Kij, a) {
  Xkij = X*Kij[,a]
  X0 = Xkij
  I = which(Xkij[,1]!=0)
  X0[I,] = standardize(Xkij[I,])$X
  return(X0)
}

EpsilonTildeEstimate = function(X, Z, p.hat, weight, ej, emax = 3000) {
  if (is.null(weight)) {
    w = weight.calculate.ps(Z, ps = p.hat, standardize = FALSE)$weight
  } else {
    w = weight
  }
  if (is.null(dim(X))) {
    X = as.matrix(X)
  } else {
    X = X
  }
  X.fit = list()
  strataIndex = list()
  for(i in 1:nrow(ej)) {
    I = which(p.hat >= ej[i,1] & p.hat < ej[i,2])
    strataIndex[[i]] = I
    if (length(I) == 1) {
      #X0 = standardize(matrix(X[I,],nrow = 1))$X
      X0 = matrix(X[I,],nrow = 1)
    } else {
      X0 = standardize(as.matrix(X[I,]))$X
    }
    X.fit[[i]] = cbind(Z[I], w[I], X0)
  }
  sample.size = sample.size.calculate(Z = Z, p.hat = p.hat, ej = ej)
  epsilon_tilde_avg = matrix(nrow = nrow(ej), ncol = ncol(X))
  sigma = matrix(nrow = nrow(ej), ncol = 2*ncol(X))
  for(i in 1:nrow(ej)) {
    n = sample.size[i,1]
    n0 = sample.size[i,3]
    n1 = sample.size[i,4]
    X0 = X.fit[[i]][which(X.fit[[i]][,1] == 0), ]
    X1 = X.fit[[i]][which(X.fit[[i]][,1] == 1), ]
    sigma0 = apply(as.matrix(X0[,3:ncol(X0)]), 2, var)
    sigma1 = apply(as.matrix(X1[,3:ncol(X1)]), 2, var)
    sigma[i,] = c(sigma0, sigma1)
    w0 = X0[,2]^2
    w1 = X1[,2]^2
    if (n <= 2 | n0 <= 2 | n1 <= 2) {
      epsilon_tilde_avg[i,] = rep(NA, ncol(X))
    } else if (n0 <= 10 | n1 <= 10) {
      epsilon_tilde_avg[i,] = rep(emax, ncol(X))
    } else if (sum(is.na(w))!= 0) {
      epsilon_tilde_avg[i,] = rep(NA, ncol(X))
    } else {
      epsilon_tilde_avg[i,] = sqrt(sum(w0)*sigma0 + sum(w1)*sigma1)
    }
  }
  return(list(epsilon_tilde_avg = epsilon_tilde_avg, sigma = sigma, sample.size = sample.size))
}

BootstrapEpsilonTildeEstimate = function(X, Z, p.hat, weight, ej, nboot, emax = 3000) {
  if (is.null(weight)) {
    #w = 1/(Z*p.hat + (1-Z)*p.hat)
    w = weight.calculate.ps(Z, ps = p.hat, standardize = FALSE)$weight
  } else {
    w = weight
  }
  X.fit = list()
  strataIndex = list()
  for(i in 1:nrow(ej)) {
    I = which(p.hat >= ej[i,1] & p.hat < ej[i,2])
    strataIndex[[i]] = I
    if (length(I) == 1) {
      X0 = standardize(matrix(X[I,],nrow = 1))$X
    } else {
      X0 = standardize(X[I,])$X
    }
    X.fit[[i]] = cbind(Z[I], w[I], X0)
  }
  sample.size = sample.size.calculate(Z = Z, p.hat = p.hat, ej = ej)
  imbalance_var = imbalance_mean = matrix(nrow = nrow(ej), ncol = ncol(X))
  imbalance_out = list()
  #sigma = matrix(nrow = nrow(ej), ncol = 2*ncol(X))
  for(i in 1:nrow(ej)) {
    n = sample.size[i,1]
    n0 = sample.size[i,3]
    n1 = sample.size[i,4]
    if (n0 <= 5 | n1 <= 5) {
      imbalance_var[i,] = imbalance_mean[i,] = rep(NA, ncol(X))
    } else if ((n0 > 5 & n0 <= 10) | (n1 >5 & n1 <= 10)) {
      imbalance_var[i,] = imbalance_mean[i,] = rep(emax, ncol(X))
    } else if (sum(is.na(w))!= 0) {
      imbalance_var[i,] = imbalance_mean[i,] = rep(NA, ncol(X))
    } else {
      X0 = X.fit[[i]][which(X.fit[[i]][,1] == 0), ]
      X1 = X.fit[[i]][which(X.fit[[i]][,1] == 1), ]
      set.seed(5050)
      imbalance = matrix(nrow = nboot, ncol = ncol(X))
      for(j in 1:nboot) {
        I0 = sample(seq(1,n0), n0, replace = TRUE)
        I1 = sample(seq(1,n1), n1, replace = TRUE)
        X0j = X0[I0,]
        X1j = X1[I1,]
        #imbalance[j,] = abs(colSums(X1j[,2]*X1j[,3:ncol(X1j)]) - colSums(X0j[,2]*X0j[,3:ncol(X0j)]))
        imbalance[j,] = colSums(X1j[,2]*X1j[,3:ncol(X1j)]) - colSums(X0j[,2]*X0j[,3:ncol(X0j)])
      }
      imbalance_var[i,] = apply(abs(imbalance), 2, var)
      imbalance_mean[i,] = apply(abs(imbalance), 2, mean)
      imbalance_out[[i]] = imbalance
    }
  }
  imbalance_sd = sqrt(imbalance_var)
  CI95.lower = imbalance_mean - 1.96*imbalance_sd
  CI95.upper = imbalance_mean + 1.96*imbalance_sd
  return(list(imbalance_var = imbalance_var, imbalance_sd = imbalance_sd,
              imbalance_mean = imbalance_mean, CI95.lower = CI95.lower,
              CI95.upper = CI95.upper, sample.size = sample.size, imbalance = imbalance_out))
}

EffectiveSampleSize = function(Z, p.hat, ej) {
  sample.size = sample.size.calculate(Z = Z, p.hat = p.hat, ej = ej)
  out = sample.size[,3]*(sample.size[,4])/(sample.size[,3]+sample.size[,4])
  return(out)
}

MaxEpsilonFind = function(X, Z, p.hat, weight = NULL, ej = ej, a) {
  sigma = EpsilonTildeEstimate(X = X, Z = Z, p.hat = p.hat, weight = weight, ej = ej)
  max.sigma = apply(sigma$epsilon_tilde_avg, 1, function(x) max(x, na.rm = T))*a
  n = EffectiveSampleSize(Z = Z, p.hat = p.hat, ej)
  out = cbind(ej, max.sigma/n)
  return(out)
}

sample.in.strata = function(p.hat, ej, In) {
  n = length(p.hat)
  strataIndex = list()
  for (i in 1:nrow(ej)) {
    I = which(p.hat >= ej[i,1] & p.hat < ej[i,2])
    strataIndex[[i]] = I
  }
  I.all = seq(1,nrow(ej))
  I.include = I.all[-In]
  if (length(I.include) == 1) {
    sample.include = strataIndex[[I.include]]
  } else if (length(I.include) != 1) {
    sample.include = strataIndex[[I.include[1]]]
    for (i in 2:length(I.include)) {
      sample.include = c(sample.include, strataIndex[[I.include[i]]])
    }
    sample.include = unique(sample.include)
    sample.include = sort(sample.include)
  }
  if (length(sample.include) == n) {
    sample.exclude = seq(1,n)[-sample.include]
  } else if (length(sample.include) == 0) {
    sample.exclude = seq(1,n)
  } else {
    sample.exclude = seq(1,n)[-sample.include]
  }
  out = list(sample.include = sample.include, no.include = length(sample.include),
             sample.exclude = sample.exclude, no.exclude = length(sample.exclude))
  return(out)
}



