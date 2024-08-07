# This file provides open source software to support the methodologies of the paper
# Branch and Bound to Assess Stability of Regression Coefficients in Uncertain Models
# by Knaeble, Hughes, Rudolph, Abramson, and Razo

# See README for details

# Helper function to calculate residuals of z
calc.z.res <- function(s, I.w, I.z) {
  z.res <- matrix(nrow = nrow(s), ncol = length(I.z))
  if (length(I.w) == 0) { # Check for root node
    # Create matrix of residuals from centering vectors of z
    return(data.frame(apply(s,2,function(x) x - mean(x))))
  } else if (length(I.z) == 0) { # Check for leaf node
    return(data.frame(matrix(nrow = nrow(s), ncol = 0)))
  } else { # Intermediate node
    # Create matrix of residuals from regressing vectors of z onto span of w
    for (z in I.z) {
      z.res[,which(I.z == z)] <- lm(s[,z]~., data = s[I.w])$residuals
    }
  }
  return(data.frame(z.res))
}

# Confounding interval function supplied by Knaeble
# Found at https://github.com/bknaeble/ConfoundingIntervals/tree/master
f=function(p,sr,lx2,ux2,ly2,uy2,lxy,uxy) {
  tol=.00001
  lx=sqrt(lx2)
  ux=sqrt(ux2)
  ly=sqrt(ly2)
  uy=sqrt(uy2)
  v=numeric(88*3)
  M=matrix(v,ncol=3)
  bx=c(lx,ux) 
  by=c(ly,uy)
  bxy=c(lxy,uxy)
  s1=function(bx,by,bxy) c(((-2*p+sqrt((2*p)^2-4*(-by*bxy)^2))/(-2*by*bxy))^2,by^2,bxy)
  s2=function(bx,by,bxy) c(((-2*p-sqrt((2*p)^2-4*(-by*bxy)^2))/(-2*by*bxy))^2,by^2,bxy)
  s3=function(bx,by,bxy) c((p+1)/(bxy+1),(p+1)/(bxy+1),bxy)
  s4=function(bx,by,bxy) c((p-1)/(bxy-1),(p-1)/(bxy-1),bxy)
  s5=function(bx,by,bxy) c(bx^2,by^2,bxy)
  s6=function(bx,by,bxy) c(bx,by,(p+sqrt(1-bx^2)*sqrt(1-by^2))/(bx*by))
  s7=function(bx,by,bxy) c(bx,by,(p-sqrt(1-bx^2)*sqrt(1-by^2))/(bx*by))
  s8=function(bx,by,bxy) c(bx^2,((-(-2*bx*bxy*p)+sqrt((-2*bx*bxy*p)^2-4*(bx^2*bxy^2+1-bx^2)*(bx^2-1+p^2)))/(2*(bx^2*bxy^2+1-bx^2)))^2,bxy)
  s9=function(bx,by,bxy) c(bx^2,((-(-2*bx*bxy*p)-sqrt((-2*bx*bxy*p)^2-4*(bx^2*bxy^2+1-bx^2)*(bx^2-1+p^2)))/(2*(bx^2*bxy^2+1-bx^2)))^2,bxy)
  s10=function(bx,by,bxy) c(((-(-2*bx*bxy*p)+sqrt((-2*bx*bxy*p)^2-4*(bx^2*bxy^2+1-bx^2)*(bx^2-1+p^2)))/(2*(bx^2*bxy^2+1-bx^2)))^2,by^2,bxy)
  s11=function(bx,by,bxy) c(((-(-2*bx*bxy*p)-sqrt((-2*bx*bxy*p)^2-4*(bx^2*bxy^2+1-bx^2)*(bx^2-1+p^2)))/(2*(bx^2*bxy^2+1-bx^2)))^2,by^2,bxy)
  w=c(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11)
  for (i in 1:2) {
    for (j in 1:2) {
      for (k in 1:2) {
        for (l in 1:11) {
          r=((l-1)*8)+(4*(i-1)+2*(j-1)+1*(k-1)+1)
          M[r,]=w[[l]](bx[i],by[j],bxy[k])
        }}}}
  M
  feas=function(arg) 
  {
    arg[1]>=(lx2-tol) & 
      arg[1]<=(ux2+tol) &
      arg[2]>=(ly2-tol) &
      arg[2]<=(uy2+tol) &
      arg[3]>=(lxy-tol) &
      arg[3]<=(uxy+tol) 
  }
  fs=apply(M,1,feas)
  N=M[which(fs==TRUE),]
  obj=function(arg) sr*(p-sqrt(arg[1])*sqrt(arg[2])*arg[3])/(1-arg[1])
  l=min(apply(N,1,obj), na.rm = TRUE)
  u=max(apply(N,1,obj), na.rm = TRUE)
  return(c(l,u))
}

# Branch and Bound Algorithm             
BB.confound.reorder <- function(x,y,s) {
  start <- Sys.time()
  # Initialize "Node" class
  setClass("Node", slots = list(beta = "ANY", x.res = "ANY", y.res = "ANY", z.res = "ANY", I.w = "ANY", I.z = "ANY"))

  # Initialize variables
  s <- data.frame(s)
  nodes.visited <- 0
  temp.min <- temp.max <- lm(y ~ x)$coef[2]

  # Initial node
  queue <- list(new("Node",
                    x.res = x - mean(x),
                    y.res = y - mean(y),
                    z.res = calc.z.res(s,numeric(0),1:ncol(s)),
                    I.w = numeric(0),
                    I.z = 1:ncol(s))
           )
  
  while (length(queue) > 0) {
    nodes.visited <- nodes.visited + 1
    # Progress update
    if (nodes.visited %% 2000 == 0) {
      print(paste0("Progress (", Sys.time(), "): ", nodes.visited, " visited. ", length(queue), " in queue."))
    }
    # Take first node out of the queue
    current <- queue[[1]]
    queue <- queue[-1]
    
    # Calculate beta.x for current node if not already done
    if (is.null(current@beta)) {
      current@beta <- as.numeric(lm(current@y.res ~ current@x.res)$coef[2])
      temp.min <- min(temp.min, current@beta)
      temp.max <- max(temp.max, current@beta)
    }

    # Find R^2 bounds
    if (length(current@I.z) > 0) {
      upper.r2.x <- summary(lm(current@x.res~., data = current@z.res))$r.squared
      upper.r2.y <- summary(lm(current@y.res~., data = current@z.res))$r.squared
    } else {
      upper.r2.x <- upper.r2.y <- 0
    }
    # Compute confounding interval
    b <- suppressWarnings(f(cor(current@x.res,current@y.res), sd(current@y.res) / sd(current@x.res), 0, upper.r2.x, 0, upper.r2.y, -1, 1))
    if (b[1] < temp.min || b[2] > temp.max) {
      # Add child nodes to queue if current is not leaf
      if (length(current@I.z) > 0) {

        # Calculate z.star
        product <- abs(cor(current@x.res, current@z.res) * cor(current@y.res, current@z.res))
        z.star <- which(product == max(product))

        # Add child nodes to queue
        queue <- c(queue,
                   new("Node",
                       beta = current@beta,
                       x.res = current@x.res,
                       y.res = current@y.res,
                       z.res = current@z.res[-z.star],
                       I.w = current@I.w,
                       I.z = current@I.z[-z.star]),
                   new("Node",
                       x.res = lm(x~., data = s[c(current@I.w, current@I.z[z.star])])$residuals,
                       y.res = lm(y~., data = s[c(current@I.w, current@I.z[z.star])])$residuals,
                       z.res = calc.z.res(s, c(current@I.w, current@I.z[z.star]), current@I.z[-z.star]),
                       I.w = c(current@I.w, current@I.z[z.star]),
                       I.z = current@I.z[-z.star]
                       )
                   )
      }
    }
  }
  end <- Sys.time()
  # Output results
  print(paste0("Nodes Visited: ", nodes.visited, " of ", (2^(ncol(s)+1)-1), " (", round(nodes.visited / (2^(ncol(s)+1)-1) * 100,2), "%)"))
  print(paste("Time elapsed: ", round(as.numeric(difftime(end, start, units = "secs")),5), "seconds"))
  return(c(temp.min, temp.max))
}
