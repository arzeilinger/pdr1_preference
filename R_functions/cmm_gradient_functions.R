#################################################################################################
#### Consumer Movement Model
#### Gradient function for Conditional Probability models NLL functions

#### Individual gradient equations as functions
# Gradient equation for p1
gr.p1 <- function(par, y){
  p1 <- par[1]
  p2 <- par[2]
  mu1 <- par[3]
  mu2 <- par[4]
  t <- y$t
  tm1 <- c(0, t[1:length(t)-1])
  tau <- t - tm1
  n1 <- y$n1
  n2 <- y$n2
  n3 <- y$n3
  n1m1 <- c(0, n1[1:length(n1)-1])
  n2m1 <- c(0, n2[1:length(n2)-1])
  N <- y$N
  gr.p1f <- ((n1*(-((4*mu2^2*p1)/(mu2*p1 + mu1*(mu2 + p2))^2) + (4*mu2)/(mu2*p1 + mu1*(mu2 + p2)) - 
                    (2*(mu1 - mu2 + p1 + p2)*(mu1^2*mu2*n1m1 - mu1*mu2^2*n1m1 - mu1*mu2*N*p1 + mu2^2*N*p1 + 
                                                2*mu1*mu2*n1m1*p1 - mu2^2*n1m1*p1 + 2*mu1*mu2*n2m1*p1 - mu2*N*p1^2 + mu2*n1m1*p1^2 + 
                                                2*mu2*n2m1*p1^2 + mu1^2*n1m1*p2 - 2*mu1*mu2*n1m1*p2 - 2*mu1*N*p1*p2 + mu2*N*p1*p2 + 
                                                mu1*n1m1*p1*p2 - mu2*n1m1*p1*p2 + 2*mu1*n2m1*p1*p2 - mu1*n1m1*p2^2 + 
                                                mu1*mu2*n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) - 
                                                mu2*N*p1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                                                mu2*n1m1*p1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                                                mu1*n1m1*p2*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))/
                    exp((1/2)*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                            4*(mu2*p1 + mu1*(mu2 + p2))))*tau)/(N*(mu2*p1 + mu1*(mu2 + p2))*
                                                                                                  ((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))^(3/2)) - 
                    (2*exp((1/2)*(-mu1 - mu2 - p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                4*(mu2*p1 + mu1*(mu2 + p2))))*tau)*(mu1 - mu2 + p1 + p2)*
                       ((-mu1^2)*mu2*n1m1 + mu1*mu2^2*n1m1 + mu1*mu2*N*p1 - mu2^2*N*p1 - 2*mu1*mu2*n1m1*p1 + 
                          mu2^2*n1m1*p1 - 2*mu1*mu2*n2m1*p1 + mu2*N*p1^2 - mu2*n1m1*p1^2 - 2*mu2*n2m1*p1^2 - 
                          mu1^2*n1m1*p2 + 2*mu1*mu2*n1m1*p2 + 2*mu1*N*p1*p2 - mu2*N*p1*p2 - mu1*n1m1*p1*p2 + 
                          mu2*n1m1*p1*p2 - 2*mu1*n2m1*p1*p2 + mu1*n1m1*p2^2 + 
                          mu1*mu2*n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) - 
                          mu2*N*p1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                          mu2*n1m1*p1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                          mu1*n1m1*p2*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))/
                    (N*(mu2*p1 + mu1*(mu2 + p2))*((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))^
                       (3/2)) + 
                    (2*(-2 + (mu2*(mu1 + mu2 + p1 + p2 - sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + 
                                                                                             mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                          (1/(N*p1))*(mu1*n1m1 - mu2*n1m1 + n1m1*p1 + 2*n2m1*p1 - n1m1*p2 + 
                                        n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))))/
                    exp((1/2)*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                            4*(mu2*p1 + mu1*(mu2 + p2))))*tau)/sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                                                      4*(mu2*p1 + mu1*(mu2 + p2))) + 
                    (2*exp((1/2)*(-mu1 - mu2 - p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                4*(mu2*p1 + mu1*(mu2 + p2))))*tau)*
                       (2 - (mu2*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                               4*(mu2*p1 + mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                          (1/(N*p1))*((-mu1)*n1m1 + mu2*n1m1 - n1m1*p1 - 2*n2m1*p1 + n1m1*p2 + 
                                        n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))))/
                    sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                    (2*exp((1/2)*(-mu1 - mu2 - p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                4*(mu2*p1 + mu1*(mu2 + p2))))*tau)*p1*
                       ((mu2*(-1 + (-mu1 + mu2 - p1 - p2)/sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + 
                                                                                              mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                          (mu2^2*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                               4*(mu2*p1 + mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2))^2 + 
                          (1/(N*p1^2))*(mu1*n1m1 - mu2*n1m1 + n1m1*p1 + 2*n2m1*p1 - n1m1*p2 - 
                                          n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))) - 
                          (2*(n1m1 + n2m1 + (1/2)*n1m1*(-1 + (-mu1 + mu2 - p1 - p2)/sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                                           4*(mu2*p1 + mu1*(mu2 + p2))))))/(N*p1)))/
                    sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                    (2*p1*(-((mu2*(-1 + (mu1 - mu2 + p1 + p2)/sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                     4*(mu2*p1 + mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2))) - 
                             (mu2^2*(mu1 + mu2 + p1 + p2 - sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + 
                                                                                               mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2))^2 - 
                             (1/(N*p1^2))*(mu1*n1m1 - mu2*n1m1 + n1m1*p1 + 2*n2m1*p1 - n1m1*p2 + 
                                             n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))) + 
                             (2*(n1m1 + n2m1 + (1/2)*n1m1*(-1 + (mu1 - mu2 + p1 + p2)/sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                                             4*(mu2*p1 + mu1*(mu2 + p2))))))/(N*p1)))/
                    exp((1/2)*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                            4*(mu2*p1 + mu1*(mu2 + p2))))*tau)/sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                                                      4*(mu2*p1 + mu1*(mu2 + p2))) + 
                    (p1*(-1 + (-mu1 + mu2 - p1 - p2)/sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                            4*(mu2*p1 + mu1*(mu2 + p2))))*
                       (-2 + (mu2*(mu1 + mu2 + p1 + p2 - sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + 
                                                                                             mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                          (1/(N*p1))*(mu1*n1m1 - mu2*n1m1 + n1m1*p1 + 2*n2m1*p1 - n1m1*p2 + 
                                        n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))*tau)/
                    exp((1/2)*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                            4*(mu2*p1 + mu1*(mu2 + p2))))*tau)/sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                                                      4*(mu2*p1 + mu1*(mu2 + p2))) + 
                    (exp((1/2)*(-mu1 - mu2 - p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                              4*(mu2*p1 + mu1*(mu2 + p2))))*tau)*p1*(-1 + (mu1 - mu2 + p1 + p2)/
                                                                                                       sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                       (2 - (mu2*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                               4*(mu2*p1 + mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                          (1/(N*p1))*((-mu1)*n1m1 + mu2*n1m1 - n1m1*p1 - 2*n2m1*p1 + n1m1*p2 + 
                                        n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))*tau)/
                    sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))/
               (2*p1*((2*mu2)/(mu2*p1 + mu1*(mu2 + p2)) + 
                        (-2 + (mu2*(mu1 + mu2 + p1 + p2 - sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                 4*(mu2*p1 + mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                           (1/(N*p1))*(mu1*n1m1 - mu2*n1m1 + n1m1*p1 + 2*n2m1*p1 - n1m1*p2 + 
                                         n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))/
                        exp((1/2)*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                4*(mu2*p1 + mu1*(mu2 + p2))))*tau)/sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                                                          4*(mu2*p1 + mu1*(mu2 + p2))) + 
                        (exp((1/2)*(-mu1 - mu2 - p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                  4*(mu2*p1 + mu1*(mu2 + p2))))*tau)*
                           (2 - (mu2*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                   4*(mu2*p1 + mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                              (1/(N*p1))*((-mu1)*n1m1 + mu2*n1m1 - n1m1*p1 - 2*n2m1*p1 + n1m1*p2 + 
                                            n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))))/
                        sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))) + 
               (n2*(-((8*mu1*mu2*p2)/(mu2*p1 + mu1*(mu2 + p2))^2) - 
                      (2*(mu1 - mu2 + p1 + p2)*(-mu1 + mu2 - p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                              4*(mu2*p1 + mu1*(mu2 + p2))))*(mu1^2*mu2*n1m1 - mu1*mu2^2*n1m1 - mu1*mu2*N*p1 + 
                                                                                                               mu2^2*N*p1 + 2*mu1*mu2*n1m1*p1 - mu2^2*n1m1*p1 + 2*mu1*mu2*n2m1*p1 - mu2*N*p1^2 + 
                                                                                                               mu2*n1m1*p1^2 + 2*mu2*n2m1*p1^2 + mu1^2*n1m1*p2 - 2*mu1*mu2*n1m1*p2 - 2*mu1*N*p1*p2 + 
                                                                                                               mu2*N*p1*p2 + mu1*n1m1*p1*p2 - mu2*n1m1*p1*p2 + 2*mu1*n2m1*p1*p2 - mu1*n1m1*p2^2 + 
                                                                                                               mu1*mu2*n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) - 
                                                                                                               mu2*N*p1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                                                                                                               mu2*n1m1*p1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                                                                                                               mu1*n1m1*p2*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))/
                      exp((1/2)*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                              4*(mu2*p1 + mu1*(mu2 + p2))))*tau)/(N*p1*(mu2*p1 + mu1*(mu2 + p2))*
                                                                                                    ((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))^(3/2)) + 
                      (2*exp((1/2)*(-mu1 - mu2 - p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                  4*(mu2*p1 + mu1*(mu2 + p2))))*tau)*(mu1 - mu2 + p1 + p2)*
                         (mu1 - mu2 + p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                         ((-mu1^2)*mu2*n1m1 + mu1*mu2^2*n1m1 + mu1*mu2*N*p1 - mu2^2*N*p1 - 2*mu1*mu2*n1m1*p1 + 
                            mu2^2*n1m1*p1 - 2*mu1*mu2*n2m1*p1 + mu2*N*p1^2 - mu2*n1m1*p1^2 - 2*mu2*n2m1*p1^2 - 
                            mu1^2*n1m1*p2 + 2*mu1*mu2*n1m1*p2 + 2*mu1*N*p1*p2 - mu2*N*p1*p2 - mu1*n1m1*p1*p2 + 
                            mu2*n1m1*p1*p2 - 2*mu1*n2m1*p1*p2 + mu1*n1m1*p2^2 + 
                            mu1*mu2*n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) - 
                            mu2*N*p1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                            mu2*n1m1*p1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                            mu1*n1m1*p2*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))/
                      (N*p1*(mu2*p1 + mu1*(mu2 + p2))*((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))^
                         (3/2)) - (2*(1 + (-mu1 + mu2 - p1 - p2)/sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                        4*(mu2*p1 + mu1*(mu2 + p2))))*
                                     (-2 + (mu2*(mu1 + mu2 + p1 + p2 - sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + 
                                                                                                           mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                                        (1/(N*p1))*(mu1*n1m1 - mu2*n1m1 + n1m1*p1 + 2*n2m1*p1 - n1m1*p2 + 
                                                      n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))))/
                      exp((1/2)*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                              4*(mu2*p1 + mu1*(mu2 + p2))))*tau)/sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                                                        4*(mu2*p1 + mu1*(mu2 + p2))) - 
                      (2*exp((1/2)*(-mu1 - mu2 - p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                  4*(mu2*p1 + mu1*(mu2 + p2))))*tau)*(1 + (mu1 - mu2 + p1 + p2)/
                                                                                                        sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                         (2 - (mu2*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                 4*(mu2*p1 + mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                            (1/(N*p1))*((-mu1)*n1m1 + mu2*n1m1 - n1m1*p1 - 2*n2m1*p1 + n1m1*p2 + 
                                          n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))))/
                      sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) - 
                      (2*exp((1/2)*(-mu1 - mu2 - p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                  4*(mu2*p1 + mu1*(mu2 + p2))))*tau)*(mu1 - mu2 + p1 - p2 + 
                                                                                                        sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                         ((mu2*(-1 + (-mu1 + mu2 - p1 - p2)/sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + 
                                                                                                mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                            (mu2^2*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                 4*(mu2*p1 + mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2))^2 + 
                            (1/(N*p1^2))*(mu1*n1m1 - mu2*n1m1 + n1m1*p1 + 2*n2m1*p1 - n1m1*p2 - 
                                            n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))) - 
                            (2*(n1m1 + n2m1 + (1/2)*n1m1*(-1 + (-mu1 + mu2 - p1 - p2)/sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                                             4*(mu2*p1 + mu1*(mu2 + p2))))))/(N*p1)))/
                      sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) - 
                      (2*(mu1 - mu2 + p1 - p2 - sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                         (-((mu2*(-1 + (mu1 - mu2 + p1 + p2)/sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                    4*(mu2*p1 + mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2))) - 
                            (mu2^2*(mu1 + mu2 + p1 + p2 - sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + 
                                                                                              mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2))^2 - 
                            (1/(N*p1^2))*(mu1*n1m1 - mu2*n1m1 + n1m1*p1 + 2*n2m1*p1 - n1m1*p2 + 
                                            n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))) + 
                            (2*(n1m1 + n2m1 + (1/2)*n1m1*(-1 + (mu1 - mu2 + p1 + p2)/sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                                            4*(mu2*p1 + mu1*(mu2 + p2))))))/(N*p1)))/
                      exp((1/2)*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                              4*(mu2*p1 + mu1*(mu2 + p2))))*tau)/sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                                                        4*(mu2*p1 + mu1*(mu2 + p2))) - 
                      ((-1 + (-mu1 + mu2 - p1 - p2)/sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                         (mu1 - mu2 + p1 - p2 - sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                         (-2 + (mu2*(mu1 + mu2 + p1 + p2 - sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + 
                                                                                               mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                            (1/(N*p1))*(mu1*n1m1 - mu2*n1m1 + n1m1*p1 + 2*n2m1*p1 - n1m1*p2 + 
                                          n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))*tau)/
                      exp((1/2)*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                              4*(mu2*p1 + mu1*(mu2 + p2))))*tau)/sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                                                        4*(mu2*p1 + mu1*(mu2 + p2))) - 
                      (exp((1/2)*(-mu1 - mu2 - p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                4*(mu2*p1 + mu1*(mu2 + p2))))*tau)*(-1 + (mu1 - mu2 + p1 + p2)/
                                                                                                      sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                         (mu1 - mu2 + p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                         (2 - (mu2*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                 4*(mu2*p1 + mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                            (1/(N*p1))*((-mu1)*n1m1 + mu2*n1m1 - n1m1*p1 - 2*n2m1*p1 + n1m1*p2 + 
                                          n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))*tau)/
                      sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))/
               (8*((mu1*p2)/(mu2*p1 + mu1*(mu2 + p2)) - 
                     ((mu1 - mu2 + p1 - p2 - sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                        (-2 + (mu2*(mu1 + mu2 + p1 + p2 - sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + 
                                                                                              mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                           (1/(N*p1))*(mu1*n1m1 - mu2*n1m1 + n1m1*p1 + 2*n2m1*p1 - n1m1*p2 + 
                                         n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))))/
                     exp((1/2)*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                             4*(mu2*p1 + mu1*(mu2 + p2))))*tau)/
                     (4*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))) - 
                     (exp((1/2)*(-mu1 - mu2 - p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                               4*(mu2*p1 + mu1*(mu2 + p2))))*tau)*(mu1 - mu2 + p1 - p2 + 
                                                                                                     sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                        (2 - (mu2*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                4*(mu2*p1 + mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                           (1/(N*p1))*((-mu1)*n1m1 + mu2*n1m1 - n1m1*p1 - 2*n2m1*p1 + n1m1*p2 + 
                                         n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))))/
                     (4*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))) + 
               (n3*((8*mu2^2*p1)/(mu2*p1 + mu1*(mu2 + p2))^2 + (8*mu1*mu2*p2)/(mu2*p1 + mu1*(mu2 + p2))^2 - 
                      (8*mu2)/(mu2*p1 + mu1*(mu2 + p2)) + 
                      (4*(mu1 - mu2 + p1 + p2)*(mu1^2*mu2*n1m1 - mu1*mu2^2*n1m1 - mu1*mu2*N*p1 + mu2^2*N*p1 + 
                                                  2*mu1*mu2*n1m1*p1 - mu2^2*n1m1*p1 + 2*mu1*mu2*n2m1*p1 - mu2*N*p1^2 + mu2*n1m1*p1^2 + 
                                                  2*mu2*n2m1*p1^2 + mu1^2*n1m1*p2 - 2*mu1*mu2*n1m1*p2 - 2*mu1*N*p1*p2 + mu2*N*p1*p2 + 
                                                  mu1*n1m1*p1*p2 - mu2*n1m1*p1*p2 + 2*mu1*n2m1*p1*p2 - mu1*n1m1*p2^2 + 
                                                  mu1*mu2*n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) - 
                                                  mu2*N*p1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                                                  mu2*n1m1*p1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                                                  mu1*n1m1*p2*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))/
                      exp((1/2)*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                              4*(mu2*p1 + mu1*(mu2 + p2))))*tau)/(N*(mu2*p1 + mu1*(mu2 + p2))*
                                                                                                    ((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))^(3/2)) + 
                      (2*(mu1 - mu2 + p1 + p2)*(-mu1 + mu2 - p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                              4*(mu2*p1 + mu1*(mu2 + p2))))*(mu1^2*mu2*n1m1 - mu1*mu2^2*n1m1 - mu1*mu2*N*p1 + 
                                                                                                               mu2^2*N*p1 + 2*mu1*mu2*n1m1*p1 - mu2^2*n1m1*p1 + 2*mu1*mu2*n2m1*p1 - mu2*N*p1^2 + 
                                                                                                               mu2*n1m1*p1^2 + 2*mu2*n2m1*p1^2 + mu1^2*n1m1*p2 - 2*mu1*mu2*n1m1*p2 - 2*mu1*N*p1*p2 + 
                                                                                                               mu2*N*p1*p2 + mu1*n1m1*p1*p2 - mu2*n1m1*p1*p2 + 2*mu1*n2m1*p1*p2 - mu1*n1m1*p2^2 + 
                                                                                                               mu1*mu2*n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) - 
                                                                                                               mu2*N*p1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                                                                                                               mu2*n1m1*p1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                                                                                                               mu1*n1m1*p2*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))/
                      exp((1/2)*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                              4*(mu2*p1 + mu1*(mu2 + p2))))*tau)/(N*p1*(mu2*p1 + mu1*(mu2 + p2))*
                                                                                                    ((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))^(3/2)) + 
                      (4*exp((1/2)*(-mu1 - mu2 - p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                  4*(mu2*p1 + mu1*(mu2 + p2))))*tau)*(mu1 - mu2 + p1 + p2)*
                         ((-mu1^2)*mu2*n1m1 + mu1*mu2^2*n1m1 + mu1*mu2*N*p1 - mu2^2*N*p1 - 2*mu1*mu2*n1m1*p1 + 
                            mu2^2*n1m1*p1 - 2*mu1*mu2*n2m1*p1 + mu2*N*p1^2 - mu2*n1m1*p1^2 - 2*mu2*n2m1*p1^2 - 
                            mu1^2*n1m1*p2 + 2*mu1*mu2*n1m1*p2 + 2*mu1*N*p1*p2 - mu2*N*p1*p2 - mu1*n1m1*p1*p2 + 
                            mu2*n1m1*p1*p2 - 2*mu1*n2m1*p1*p2 + mu1*n1m1*p2^2 + 
                            mu1*mu2*n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) - 
                            mu2*N*p1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                            mu2*n1m1*p1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                            mu1*n1m1*p2*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))/
                      (N*(mu2*p1 + mu1*(mu2 + p2))*((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))^
                         (3/2)) - (2*exp((1/2)*(-mu1 - mu2 - p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                              4*(mu2*p1 + mu1*(mu2 + p2))))*tau)*(mu1 - mu2 + p1 + p2)*
                                     (mu1 - mu2 + p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                                     ((-mu1^2)*mu2*n1m1 + mu1*mu2^2*n1m1 + mu1*mu2*N*p1 - mu2^2*N*p1 - 2*mu1*mu2*n1m1*p1 + 
                                        mu2^2*n1m1*p1 - 2*mu1*mu2*n2m1*p1 + mu2*N*p1^2 - mu2*n1m1*p1^2 - 2*mu2*n2m1*p1^2 - 
                                        mu1^2*n1m1*p2 + 2*mu1*mu2*n1m1*p2 + 2*mu1*N*p1*p2 - mu2*N*p1*p2 - mu1*n1m1*p1*p2 + 
                                        mu2*n1m1*p1*p2 - 2*mu1*n2m1*p1*p2 + mu1*n1m1*p2^2 + 
                                        mu1*mu2*n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) - 
                                        mu2*N*p1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                                        mu2*n1m1*p1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                                        mu1*n1m1*p2*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))/
                      (N*p1*(mu2*p1 + mu1*(mu2 + p2))*((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))^
                         (3/2)) - 
                      (4*(-2 + (mu2*(mu1 + mu2 + p1 + p2 - sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + 
                                                                                               mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                            (1/(N*p1))*(mu1*n1m1 - mu2*n1m1 + n1m1*p1 + 2*n2m1*p1 - n1m1*p2 + 
                                          n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))))/
                      exp((1/2)*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                              4*(mu2*p1 + mu1*(mu2 + p2))))*tau)/sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                                                        4*(mu2*p1 + mu1*(mu2 + p2))) + 
                      (2*(1 + (-mu1 + mu2 - p1 - p2)/sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                         (-2 + (mu2*(mu1 + mu2 + p1 + p2 - sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + 
                                                                                               mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                            (1/(N*p1))*(mu1*n1m1 - mu2*n1m1 + n1m1*p1 + 2*n2m1*p1 - n1m1*p2 + 
                                          n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))))/
                      exp((1/2)*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                              4*(mu2*p1 + mu1*(mu2 + p2))))*tau)/sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                                                        4*(mu2*p1 + mu1*(mu2 + p2))) - 
                      (4*exp((1/2)*(-mu1 - mu2 - p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                  4*(mu2*p1 + mu1*(mu2 + p2))))*tau)*
                         (2 - (mu2*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                 4*(mu2*p1 + mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                            (1/(N*p1))*((-mu1)*n1m1 + mu2*n1m1 - n1m1*p1 - 2*n2m1*p1 + n1m1*p2 + 
                                          n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))))/
                      sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                      (2*exp((1/2)*(-mu1 - mu2 - p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                  4*(mu2*p1 + mu1*(mu2 + p2))))*tau)*(1 + (mu1 - mu2 + p1 + p2)/
                                                                                                        sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                         (2 - (mu2*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                 4*(mu2*p1 + mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                            (1/(N*p1))*((-mu1)*n1m1 + mu2*n1m1 - n1m1*p1 - 2*n2m1*p1 + n1m1*p2 + 
                                          n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))))/
                      sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) - 
                      (4*exp((1/2)*(-mu1 - mu2 - p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                  4*(mu2*p1 + mu1*(mu2 + p2))))*tau)*p1*
                         ((mu2*(-1 + (-mu1 + mu2 - p1 - p2)/sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + 
                                                                                                mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                            (mu2^2*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                 4*(mu2*p1 + mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2))^2 + 
                            (1/(N*p1^2))*(mu1*n1m1 - mu2*n1m1 + n1m1*p1 + 2*n2m1*p1 - n1m1*p2 - 
                                            n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))) - 
                            (2*(n1m1 + n2m1 + (1/2)*n1m1*(-1 + (-mu1 + mu2 - p1 - p2)/sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                                             4*(mu2*p1 + mu1*(mu2 + p2))))))/(N*p1)))/
                      sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                      (2*exp((1/2)*(-mu1 - mu2 - p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                  4*(mu2*p1 + mu1*(mu2 + p2))))*tau)*(mu1 - mu2 + p1 - p2 + 
                                                                                                        sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                         ((mu2*(-1 + (-mu1 + mu2 - p1 - p2)/sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + 
                                                                                                mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                            (mu2^2*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                 4*(mu2*p1 + mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2))^2 + 
                            (1/(N*p1^2))*(mu1*n1m1 - mu2*n1m1 + n1m1*p1 + 2*n2m1*p1 - n1m1*p2 - 
                                            n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))) - 
                            (2*(n1m1 + n2m1 + (1/2)*n1m1*(-1 + (-mu1 + mu2 - p1 - p2)/sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                                             4*(mu2*p1 + mu1*(mu2 + p2))))))/(N*p1)))/
                      sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) - 
                      (4*p1*(-((mu2*(-1 + (mu1 - mu2 + p1 + p2)/sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                       4*(mu2*p1 + mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2))) - 
                               (mu2^2*(mu1 + mu2 + p1 + p2 - sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + 
                                                                                                 mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2))^2 - 
                               (1/(N*p1^2))*(mu1*n1m1 - mu2*n1m1 + n1m1*p1 + 2*n2m1*p1 - n1m1*p2 + 
                                               n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))) + 
                               (2*(n1m1 + n2m1 + (1/2)*n1m1*(-1 + (mu1 - mu2 + p1 + p2)/sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                                               4*(mu2*p1 + mu1*(mu2 + p2))))))/(N*p1)))/
                      exp((1/2)*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                              4*(mu2*p1 + mu1*(mu2 + p2))))*tau)/sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                                                        4*(mu2*p1 + mu1*(mu2 + p2))) + 
                      (2*(mu1 - mu2 + p1 - p2 - sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                         (-((mu2*(-1 + (mu1 - mu2 + p1 + p2)/sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                    4*(mu2*p1 + mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2))) - 
                            (mu2^2*(mu1 + mu2 + p1 + p2 - sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + 
                                                                                              mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2))^2 - 
                            (1/(N*p1^2))*(mu1*n1m1 - mu2*n1m1 + n1m1*p1 + 2*n2m1*p1 - n1m1*p2 + 
                                            n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))) + 
                            (2*(n1m1 + n2m1 + (1/2)*n1m1*(-1 + (mu1 - mu2 + p1 + p2)/sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                                            4*(mu2*p1 + mu1*(mu2 + p2))))))/(N*p1)))/
                      exp((1/2)*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                              4*(mu2*p1 + mu1*(mu2 + p2))))*tau)/sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                                                        4*(mu2*p1 + mu1*(mu2 + p2))) - 
                      (2*p1*(-1 + (-mu1 + mu2 - p1 - p2)/sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                4*(mu2*p1 + mu1*(mu2 + p2))))*
                         (-2 + (mu2*(mu1 + mu2 + p1 + p2 - sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + 
                                                                                               mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                            (1/(N*p1))*(mu1*n1m1 - mu2*n1m1 + n1m1*p1 + 2*n2m1*p1 - n1m1*p2 + 
                                          n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))*tau)/
                      exp((1/2)*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                              4*(mu2*p1 + mu1*(mu2 + p2))))*tau)/sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                                                        4*(mu2*p1 + mu1*(mu2 + p2))) + 
                      ((-1 + (-mu1 + mu2 - p1 - p2)/sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                         (mu1 - mu2 + p1 - p2 - sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                         (-2 + (mu2*(mu1 + mu2 + p1 + p2 - sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + 
                                                                                               mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                            (1/(N*p1))*(mu1*n1m1 - mu2*n1m1 + n1m1*p1 + 2*n2m1*p1 - n1m1*p2 + 
                                          n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))*tau)/
                      exp((1/2)*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                              4*(mu2*p1 + mu1*(mu2 + p2))))*tau)/sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                                                        4*(mu2*p1 + mu1*(mu2 + p2))) - 
                      (2*exp((1/2)*(-mu1 - mu2 - p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                  4*(mu2*p1 + mu1*(mu2 + p2))))*tau)*p1*(-1 + (mu1 - mu2 + p1 + p2)/
                                                                                                           sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                         (2 - (mu2*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                 4*(mu2*p1 + mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                            (1/(N*p1))*((-mu1)*n1m1 + mu2*n1m1 - n1m1*p1 - 2*n2m1*p1 + n1m1*p2 + 
                                          n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))*tau)/
                      sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                      (exp((1/2)*(-mu1 - mu2 - p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                4*(mu2*p1 + mu1*(mu2 + p2))))*tau)*(-1 + (mu1 - mu2 + p1 + p2)/
                                                                                                      sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                         (mu1 - mu2 + p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                         (2 - (mu2*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                 4*(mu2*p1 + mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                            (1/(N*p1))*((-mu1)*n1m1 + mu2*n1m1 - n1m1*p1 - 2*n2m1*p1 + n1m1*p2 + 
                                          n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))*tau)/
                      sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))/
               (2*(4 - (4*mu2*p1)/(mu2*p1 + mu1*(mu2 + p2)) - (4*mu1*p2)/(mu2*p1 + mu1*(mu2 + p2)) - 
                     (2*p1*(-2 + (mu2*(mu1 + mu2 + p1 + p2 - sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                    4*(mu2*p1 + mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                              (1/(N*p1))*(mu1*n1m1 - mu2*n1m1 + n1m1*p1 + 2*n2m1*p1 - n1m1*p2 + 
                                            n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))))/
                     exp((1/2)*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                             4*(mu2*p1 + mu1*(mu2 + p2))))*tau)/sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                                                       4*(mu2*p1 + mu1*(mu2 + p2))) + 
                     ((mu1 - mu2 + p1 - p2 - sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                        (-2 + (mu2*(mu1 + mu2 + p1 + p2 - sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + 
                                                                                              mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                           (1/(N*p1))*(mu1*n1m1 - mu2*n1m1 + n1m1*p1 + 2*n2m1*p1 - n1m1*p2 + 
                                         n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))))/
                     exp((1/2)*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                             4*(mu2*p1 + mu1*(mu2 + p2))))*tau)/sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                                                       4*(mu2*p1 + mu1*(mu2 + p2))) - 
                     (2*exp((1/2)*(-mu1 - mu2 - p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                 4*(mu2*p1 + mu1*(mu2 + p2))))*tau)*p1*
                        (2 - (mu2*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                4*(mu2*p1 + mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                           (1/(N*p1))*((-mu1)*n1m1 + mu2*n1m1 - n1m1*p1 - 2*n2m1*p1 + n1m1*p2 + 
                                         n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))))/
                     sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                     (exp((1/2)*(-mu1 - mu2 - p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                               4*(mu2*p1 + mu1*(mu2 + p2))))*tau)*(mu1 - mu2 + p1 - p2 + 
                                                                                                     sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                        (2 - (mu2*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                4*(mu2*p1 + mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                           (1/(N*p1))*((-mu1)*n1m1 + mu2*n1m1 - n1m1*p1 - 2*n2m1*p1 + n1m1*p2 + 
                                         n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))))/
                     sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))))
  return(gr.p1f)
}

# Gradient equation for p2
gr.p2 <- function(par, y){
  p1 <- par[1]
  p2 <- par[2]
  mu1 <- par[3]
  mu2 <- par[4]
  t <- y$t
  tm1 <- c(0, t[1:length(t)-1])
  tau <- t - tm1
  n1 <- y$n1
  n2 <- y$n2
  n3 <- y$n3
  n1m1 <- c(0, n1[1:length(n1)-1])
  n2m1 <- c(0, n2[1:length(n2)-1])
  N <- y$N
  gr.p2f <- ((n1*(-((2*mu1*mu2)/(mu2*p1 + mu1*(mu2 + p2))^2) + 
                    ((mu1 - mu2 - p1 - p2)*(mu1^2*mu2*n1m1 - mu1*mu2^2*n1m1 - mu1*mu2*N*p1 + mu2^2*N*p1 + 
                                              2*mu1*mu2*n1m1*p1 - mu2^2*n1m1*p1 + 2*mu1*mu2*n2m1*p1 - mu2*N*p1^2 + mu2*n1m1*p1^2 + 
                                              2*mu2*n2m1*p1^2 + mu1^2*n1m1*p2 - 2*mu1*mu2*n1m1*p2 - 2*mu1*N*p1*p2 + mu2*N*p1*p2 + 
                                              mu1*n1m1*p1*p2 - mu2*n1m1*p1*p2 + 2*mu1*n2m1*p1*p2 - mu1*n1m1*p2^2 + 
                                              mu1*mu2*n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) - 
                                              mu2*N*p1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                                              mu2*n1m1*p1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                                              mu1*n1m1*p2*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))/
                    exp((1/2)*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                            4*(mu2*p1 + mu1*(mu2 + p2))))*tau)/(N*p1*(mu2*p1 + mu1*(mu2 + p2))*
                                                                                                  ((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))^(3/2)) + 
                    (exp((1/2)*(-mu1 - mu2 - p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                              4*(mu2*p1 + mu1*(mu2 + p2))))*tau)*(mu1 - mu2 - p1 - p2)*
                       ((-mu1^2)*mu2*n1m1 + mu1*mu2^2*n1m1 + mu1*mu2*N*p1 - mu2^2*N*p1 - 2*mu1*mu2*n1m1*p1 + 
                          mu2^2*n1m1*p1 - 2*mu1*mu2*n2m1*p1 + mu2*N*p1^2 - mu2*n1m1*p1^2 - 2*mu2*n2m1*p1^2 - 
                          mu1^2*n1m1*p2 + 2*mu1*mu2*n1m1*p2 + 2*mu1*N*p1*p2 - mu2*N*p1*p2 - mu1*n1m1*p1*p2 + 
                          mu2*n1m1*p1*p2 - 2*mu1*n2m1*p1*p2 + mu1*n1m1*p2^2 + 
                          mu1*mu2*n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) - 
                          mu2*N*p1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                          mu2*n1m1*p1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                          mu1*n1m1*p2*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))/
                    (N*p1*(mu2*p1 + mu1*(mu2 + p2))*((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))^
                       (3/2)) + 
                    ((n1m1*(-1 + (-mu1 + mu2 + p1 + p2)/sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + 
                                                                                            mu1*(mu2 + p2)))))/(N*p1) - 
                       (mu2*(-1 + (-mu1 + mu2 + p1 + p2)/sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + 
                                                                                             mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) - 
                       (mu1*mu2*(mu1 + mu2 + p1 + p2 - sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                              4*(mu2*p1 + mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2))^2)/
                    exp((1/2)*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                            4*(mu2*p1 + mu1*(mu2 + p2))))*tau)/sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                                                      4*(mu2*p1 + mu1*(mu2 + p2))) + 
                    (exp((1/2)*(-mu1 - mu2 - p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                              4*(mu2*p1 + mu1*(mu2 + p2))))*tau)*
                       (-((n1m1*(-1 + (mu1 - mu2 - p1 - p2)/sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + 
                                                                                                mu1*(mu2 + p2)))))/(N*p1)) + 
                          (mu2*(-1 + (mu1 - mu2 - p1 - p2)/sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + 
                                                                                               mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                          (mu1*mu2*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                 4*(mu2*p1 + mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2))^2))/
                    sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                    ((-1 + (mu1 - mu2 - p1 - p2)/sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                       (-2 + (mu2*(mu1 + mu2 + p1 + p2 - sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + 
                                                                                             mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                          (1/(N*p1))*(mu1*n1m1 - mu2*n1m1 + n1m1*p1 + 2*n2m1*p1 - n1m1*p2 + 
                                        n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))*tau)/
                    exp((1/2)*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                            4*(mu2*p1 + mu1*(mu2 + p2))))*tau)/
                    (2*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))) + 
                    (exp((1/2)*(-mu1 - mu2 - p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                              4*(mu2*p1 + mu1*(mu2 + p2))))*tau)*(-1 + (-mu1 + mu2 + p1 + p2)/
                                                                                                    sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                       (2 - (mu2*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                               4*(mu2*p1 + mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                          (1/(N*p1))*((-mu1)*n1m1 + mu2*n1m1 - n1m1*p1 - 2*n2m1*p1 + n1m1*p2 + 
                                        n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))*tau)/
                    (2*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))))/
               ((2*mu2)/(mu2*p1 + mu1*(mu2 + p2)) + 
                  (-2 + (mu2*(mu1 + mu2 + p1 + p2 - sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                           4*(mu2*p1 + mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                     (1/(N*p1))*(mu1*n1m1 - mu2*n1m1 + n1m1*p1 + 2*n2m1*p1 - n1m1*p2 + 
                                   n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))/
                  exp((1/2)*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                          4*(mu2*p1 + mu1*(mu2 + p2))))*tau)/sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                                                    4*(mu2*p1 + mu1*(mu2 + p2))) + 
                  (exp((1/2)*(-mu1 - mu2 - p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                            4*(mu2*p1 + mu1*(mu2 + p2))))*tau)*
                     (2 - (mu2*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                             4*(mu2*p1 + mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                        (1/(N*p1))*((-mu1)*n1m1 + mu2*n1m1 - n1m1*p1 - 2*n2m1*p1 + n1m1*p2 + 
                                      n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))))/
                  sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))) + 
               (n2*(-((8*mu1^2*p2)/(mu2*p1 + mu1*(mu2 + p2))^2) + (8*mu1)/(mu2*p1 + mu1*(mu2 + p2)) + 
                      (2*(mu1 - mu2 - p1 - p2)*(-mu1 + mu2 - p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                              4*(mu2*p1 + mu1*(mu2 + p2))))*(mu1^2*mu2*n1m1 - mu1*mu2^2*n1m1 - mu1*mu2*N*p1 + 
                                                                                                               mu2^2*N*p1 + 2*mu1*mu2*n1m1*p1 - mu2^2*n1m1*p1 + 2*mu1*mu2*n2m1*p1 - mu2*N*p1^2 + 
                                                                                                               mu2*n1m1*p1^2 + 2*mu2*n2m1*p1^2 + mu1^2*n1m1*p2 - 2*mu1*mu2*n1m1*p2 - 2*mu1*N*p1*p2 + 
                                                                                                               mu2*N*p1*p2 + mu1*n1m1*p1*p2 - mu2*n1m1*p1*p2 + 2*mu1*n2m1*p1*p2 - mu1*n1m1*p2^2 + 
                                                                                                               mu1*mu2*n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) - 
                                                                                                               mu2*N*p1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                                                                                                               mu2*n1m1*p1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                                                                                                               mu1*n1m1*p2*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))/
                      exp((1/2)*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                              4*(mu2*p1 + mu1*(mu2 + p2))))*tau)/(N*p1*(mu2*p1 + mu1*(mu2 + p2))*
                                                                                                    ((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))^(3/2)) - 
                      (2*exp((1/2)*(-mu1 - mu2 - p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                  4*(mu2*p1 + mu1*(mu2 + p2))))*tau)*(mu1 - mu2 - p1 - p2)*
                         (mu1 - mu2 + p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                         ((-mu1^2)*mu2*n1m1 + mu1*mu2^2*n1m1 + mu1*mu2*N*p1 - mu2^2*N*p1 - 2*mu1*mu2*n1m1*p1 + 
                            mu2^2*n1m1*p1 - 2*mu1*mu2*n2m1*p1 + mu2*N*p1^2 - mu2*n1m1*p1^2 - 2*mu2*n2m1*p1^2 - 
                            mu1^2*n1m1*p2 + 2*mu1*mu2*n1m1*p2 + 2*mu1*N*p1*p2 - mu2*N*p1*p2 - mu1*n1m1*p1*p2 + 
                            mu2*n1m1*p1*p2 - 2*mu1*n2m1*p1*p2 + mu1*n1m1*p2^2 + 
                            mu1*mu2*n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) - 
                            mu2*N*p1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                            mu2*n1m1*p1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                            mu1*n1m1*p2*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))/
                      (N*p1*(mu2*p1 + mu1*(mu2 + p2))*((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))^
                         (3/2)) - (2*(mu1 - mu2 + p1 - p2 - sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                   4*(mu2*p1 + mu1*(mu2 + p2))))*
                                     ((n1m1*(-1 + (-mu1 + mu2 + p1 + p2)/sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + 
                                                                                                             mu1*(mu2 + p2)))))/(N*p1) - 
                                        (mu2*(-1 + (-mu1 + mu2 + p1 + p2)/sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + 
                                                                                                              mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) - 
                                        (mu1*mu2*(mu1 + mu2 + p1 + p2 - sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + 
                                                                                                            mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2))^2))/
                      exp((1/2)*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                              4*(mu2*p1 + mu1*(mu2 + p2))))*tau)/sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                                                        4*(mu2*p1 + mu1*(mu2 + p2))) - 
                      (2*exp((1/2)*(-mu1 - mu2 - p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                  4*(mu2*p1 + mu1*(mu2 + p2))))*tau)*(mu1 - mu2 + p1 - p2 + 
                                                                                                        sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                         (-((n1m1*(-1 + (mu1 - mu2 - p1 - p2)/sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + 
                                                                                                  mu1*(mu2 + p2)))))/(N*p1)) + 
                            (mu2*(-1 + (mu1 - mu2 - p1 - p2)/sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + 
                                                                                                 mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                            (mu1*mu2*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                   4*(mu2*p1 + mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2))^2))/
                      sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) - 
                      (2*(-1 + (mu1 - mu2 - p1 - p2)/sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                         (-2 + (mu2*(mu1 + mu2 + p1 + p2 - sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + 
                                                                                               mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                            (1/(N*p1))*(mu1*n1m1 - mu2*n1m1 + n1m1*p1 + 2*n2m1*p1 - n1m1*p2 + 
                                          n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))))/
                      exp((1/2)*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                              4*(mu2*p1 + mu1*(mu2 + p2))))*tau)/sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                                                        4*(mu2*p1 + mu1*(mu2 + p2))) - 
                      (2*exp((1/2)*(-mu1 - mu2 - p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                  4*(mu2*p1 + mu1*(mu2 + p2))))*tau)*(-1 + (-mu1 + mu2 + p1 + p2)/
                                                                                                        sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                         (2 - (mu2*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                 4*(mu2*p1 + mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                            (1/(N*p1))*((-mu1)*n1m1 + mu2*n1m1 - n1m1*p1 - 2*n2m1*p1 + n1m1*p2 + 
                                          n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))))/
                      sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) - 
                      ((-1 + (mu1 - mu2 - p1 - p2)/sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                         (mu1 - mu2 + p1 - p2 - sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                         (-2 + (mu2*(mu1 + mu2 + p1 + p2 - sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + 
                                                                                               mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                            (1/(N*p1))*(mu1*n1m1 - mu2*n1m1 + n1m1*p1 + 2*n2m1*p1 - n1m1*p2 + 
                                          n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))*tau)/
                      exp((1/2)*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                              4*(mu2*p1 + mu1*(mu2 + p2))))*tau)/sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                                                        4*(mu2*p1 + mu1*(mu2 + p2))) - 
                      (exp((1/2)*(-mu1 - mu2 - p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                4*(mu2*p1 + mu1*(mu2 + p2))))*tau)*(-1 + (-mu1 + mu2 + p1 + p2)/
                                                                                                      sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                         (mu1 - mu2 + p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                         (2 - (mu2*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                 4*(mu2*p1 + mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                            (1/(N*p1))*((-mu1)*n1m1 + mu2*n1m1 - n1m1*p1 - 2*n2m1*p1 + n1m1*p2 + 
                                          n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))*tau)/
                      sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))/
               (8*((mu1*p2)/(mu2*p1 + mu1*(mu2 + p2)) - 
                     ((mu1 - mu2 + p1 - p2 - sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                        (-2 + (mu2*(mu1 + mu2 + p1 + p2 - sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + 
                                                                                              mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                           (1/(N*p1))*(mu1*n1m1 - mu2*n1m1 + n1m1*p1 + 2*n2m1*p1 - n1m1*p2 + 
                                         n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))))/
                     exp((1/2)*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                             4*(mu2*p1 + mu1*(mu2 + p2))))*tau)/
                     (4*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))) - 
                     (exp((1/2)*(-mu1 - mu2 - p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                               4*(mu2*p1 + mu1*(mu2 + p2))))*tau)*(mu1 - mu2 + p1 - p2 + 
                                                                                                     sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                        (2 - (mu2*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                4*(mu2*p1 + mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                           (1/(N*p1))*((-mu1)*n1m1 + mu2*n1m1 - n1m1*p1 - 2*n2m1*p1 + n1m1*p2 + 
                                         n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))))/
                     (4*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))) + 
               (n3*((8*mu1*mu2*p1)/(mu2*p1 + mu1*(mu2 + p2))^2 + (8*mu1^2*p2)/(mu2*p1 + mu1*(mu2 + p2))^2 - 
                      (8*mu1)/(mu2*p1 + mu1*(mu2 + p2)) - 
                      (4*(mu1 - mu2 - p1 - p2)*(mu1^2*mu2*n1m1 - mu1*mu2^2*n1m1 - mu1*mu2*N*p1 + mu2^2*N*p1 + 
                                                  2*mu1*mu2*n1m1*p1 - mu2^2*n1m1*p1 + 2*mu1*mu2*n2m1*p1 - mu2*N*p1^2 + mu2*n1m1*p1^2 + 
                                                  2*mu2*n2m1*p1^2 + mu1^2*n1m1*p2 - 2*mu1*mu2*n1m1*p2 - 2*mu1*N*p1*p2 + mu2*N*p1*p2 + 
                                                  mu1*n1m1*p1*p2 - mu2*n1m1*p1*p2 + 2*mu1*n2m1*p1*p2 - mu1*n1m1*p2^2 + 
                                                  mu1*mu2*n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) - 
                                                  mu2*N*p1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                                                  mu2*n1m1*p1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                                                  mu1*n1m1*p2*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))/
                      exp((1/2)*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                              4*(mu2*p1 + mu1*(mu2 + p2))))*tau)/(N*(mu2*p1 + mu1*(mu2 + p2))*
                                                                                                    ((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))^(3/2)) - 
                      (2*(mu1 - mu2 - p1 - p2)*(-mu1 + mu2 - p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                              4*(mu2*p1 + mu1*(mu2 + p2))))*(mu1^2*mu2*n1m1 - mu1*mu2^2*n1m1 - mu1*mu2*N*p1 + 
                                                                                                               mu2^2*N*p1 + 2*mu1*mu2*n1m1*p1 - mu2^2*n1m1*p1 + 2*mu1*mu2*n2m1*p1 - mu2*N*p1^2 + 
                                                                                                               mu2*n1m1*p1^2 + 2*mu2*n2m1*p1^2 + mu1^2*n1m1*p2 - 2*mu1*mu2*n1m1*p2 - 2*mu1*N*p1*p2 + 
                                                                                                               mu2*N*p1*p2 + mu1*n1m1*p1*p2 - mu2*n1m1*p1*p2 + 2*mu1*n2m1*p1*p2 - mu1*n1m1*p2^2 + 
                                                                                                               mu1*mu2*n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) - 
                                                                                                               mu2*N*p1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                                                                                                               mu2*n1m1*p1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                                                                                                               mu1*n1m1*p2*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))/
                      exp((1/2)*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                              4*(mu2*p1 + mu1*(mu2 + p2))))*tau)/(N*p1*(mu2*p1 + mu1*(mu2 + p2))*
                                                                                                    ((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))^(3/2)) + 
                      (4*exp((1/2)*(-mu1 - mu2 - p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                  4*(mu2*p1 + mu1*(mu2 + p2))))*tau)*(-mu1 + mu2 + p1 + p2)*
                         ((-mu1^2)*mu2*n1m1 + mu1*mu2^2*n1m1 + mu1*mu2*N*p1 - mu2^2*N*p1 - 2*mu1*mu2*n1m1*p1 + 
                            mu2^2*n1m1*p1 - 2*mu1*mu2*n2m1*p1 + mu2*N*p1^2 - mu2*n1m1*p1^2 - 2*mu2*n2m1*p1^2 - 
                            mu1^2*n1m1*p2 + 2*mu1*mu2*n1m1*p2 + 2*mu1*N*p1*p2 - mu2*N*p1*p2 - mu1*n1m1*p1*p2 + 
                            mu2*n1m1*p1*p2 - 2*mu1*n2m1*p1*p2 + mu1*n1m1*p2^2 + 
                            mu1*mu2*n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) - 
                            mu2*N*p1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                            mu2*n1m1*p1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                            mu1*n1m1*p2*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))/
                      (N*(mu2*p1 + mu1*(mu2 + p2))*((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))^
                         (3/2)) + (2*exp((1/2)*(-mu1 - mu2 - p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                              4*(mu2*p1 + mu1*(mu2 + p2))))*tau)*(mu1 - mu2 - p1 - p2)*
                                     (mu1 - mu2 + p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                                     ((-mu1^2)*mu2*n1m1 + mu1*mu2^2*n1m1 + mu1*mu2*N*p1 - mu2^2*N*p1 - 2*mu1*mu2*n1m1*p1 + 
                                        mu2^2*n1m1*p1 - 2*mu1*mu2*n2m1*p1 + mu2*N*p1^2 - mu2*n1m1*p1^2 - 2*mu2*n2m1*p1^2 - 
                                        mu1^2*n1m1*p2 + 2*mu1*mu2*n1m1*p2 + 2*mu1*N*p1*p2 - mu2*N*p1*p2 - mu1*n1m1*p1*p2 + 
                                        mu2*n1m1*p1*p2 - 2*mu1*n2m1*p1*p2 + mu1*n1m1*p2^2 + 
                                        mu1*mu2*n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) - 
                                        mu2*N*p1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                                        mu2*n1m1*p1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                                        mu1*n1m1*p2*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))/
                      (N*p1*(mu2*p1 + mu1*(mu2 + p2))*((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))^
                         (3/2)) - 
                      (4*p1*((n1m1*(-1 + (-mu1 + mu2 + p1 + p2)/sqrt((mu1 + mu2 + p1 + p2)^2 - 4*
                                                                       (mu2*p1 + mu1*(mu2 + p2)))))/(N*p1) - 
                               (mu2*(-1 + (-mu1 + mu2 + p1 + p2)/sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + 
                                                                                                     mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) - 
                               (mu1*mu2*(mu1 + mu2 + p1 + p2 - sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + 
                                                                                                   mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2))^2))/
                      exp((1/2)*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                              4*(mu2*p1 + mu1*(mu2 + p2))))*tau)/sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                                                        4*(mu2*p1 + mu1*(mu2 + p2))) + 
                      (2*(mu1 - mu2 + p1 - p2 - sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                         ((n1m1*(-1 + (-mu1 + mu2 + p1 + p2)/sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + 
                                                                                                 mu1*(mu2 + p2)))))/(N*p1) - 
                            (mu2*(-1 + (-mu1 + mu2 + p1 + p2)/sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + 
                                                                                                  mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) - 
                            (mu1*mu2*(mu1 + mu2 + p1 + p2 - sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + 
                                                                                                mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2))^2))/
                      exp((1/2)*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                              4*(mu2*p1 + mu1*(mu2 + p2))))*tau)/sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                                                        4*(mu2*p1 + mu1*(mu2 + p2))) - 
                      (4*exp((1/2)*(-mu1 - mu2 - p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                  4*(mu2*p1 + mu1*(mu2 + p2))))*tau)*p1*
                         (-((n1m1*(-1 + (mu1 - mu2 - p1 - p2)/sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + 
                                                                                                  mu1*(mu2 + p2)))))/(N*p1)) + 
                            (mu2*(-1 + (mu1 - mu2 - p1 - p2)/sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + 
                                                                                                 mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                            (mu1*mu2*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                   4*(mu2*p1 + mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2))^2))/
                      sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                      (2*exp((1/2)*(-mu1 - mu2 - p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                  4*(mu2*p1 + mu1*(mu2 + p2))))*tau)*(mu1 - mu2 + p1 - p2 + 
                                                                                                        sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                         (-((n1m1*(-1 + (mu1 - mu2 - p1 - p2)/sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + 
                                                                                                  mu1*(mu2 + p2)))))/(N*p1)) + 
                            (mu2*(-1 + (mu1 - mu2 - p1 - p2)/sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + 
                                                                                                 mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                            (mu1*mu2*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                   4*(mu2*p1 + mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2))^2))/
                      sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                      (2*(-1 + (mu1 - mu2 - p1 - p2)/sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                         (-2 + (mu2*(mu1 + mu2 + p1 + p2 - sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + 
                                                                                               mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                            (1/(N*p1))*(mu1*n1m1 - mu2*n1m1 + n1m1*p1 + 2*n2m1*p1 - n1m1*p2 + 
                                          n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))))/
                      exp((1/2)*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                              4*(mu2*p1 + mu1*(mu2 + p2))))*tau)/sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                                                        4*(mu2*p1 + mu1*(mu2 + p2))) + 
                      (2*exp((1/2)*(-mu1 - mu2 - p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                  4*(mu2*p1 + mu1*(mu2 + p2))))*tau)*(-1 + (-mu1 + mu2 + p1 + p2)/
                                                                                                        sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                         (2 - (mu2*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                 4*(mu2*p1 + mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                            (1/(N*p1))*((-mu1)*n1m1 + mu2*n1m1 - n1m1*p1 - 2*n2m1*p1 + n1m1*p2 + 
                                          n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))))/
                      sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) - 
                      (2*p1*(-1 + (mu1 - mu2 - p1 - p2)/sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                               4*(mu2*p1 + mu1*(mu2 + p2))))*
                         (-2 + (mu2*(mu1 + mu2 + p1 + p2 - sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + 
                                                                                               mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                            (1/(N*p1))*(mu1*n1m1 - mu2*n1m1 + n1m1*p1 + 2*n2m1*p1 - n1m1*p2 + 
                                          n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))*tau)/
                      exp((1/2)*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                              4*(mu2*p1 + mu1*(mu2 + p2))))*tau)/sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                                                        4*(mu2*p1 + mu1*(mu2 + p2))) + 
                      ((-1 + (mu1 - mu2 - p1 - p2)/sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                         (mu1 - mu2 + p1 - p2 - sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                         (-2 + (mu2*(mu1 + mu2 + p1 + p2 - sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + 
                                                                                               mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                            (1/(N*p1))*(mu1*n1m1 - mu2*n1m1 + n1m1*p1 + 2*n2m1*p1 - n1m1*p2 + 
                                          n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))*tau)/
                      exp((1/2)*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                              4*(mu2*p1 + mu1*(mu2 + p2))))*tau)/sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                                                        4*(mu2*p1 + mu1*(mu2 + p2))) - 
                      (2*exp((1/2)*(-mu1 - mu2 - p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                  4*(mu2*p1 + mu1*(mu2 + p2))))*tau)*p1*(-1 + (-mu1 + mu2 + p1 + p2)/
                                                                                                           sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                         (2 - (mu2*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                 4*(mu2*p1 + mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                            (1/(N*p1))*((-mu1)*n1m1 + mu2*n1m1 - n1m1*p1 - 2*n2m1*p1 + n1m1*p2 + 
                                          n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))*tau)/
                      sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                      (exp((1/2)*(-mu1 - mu2 - p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                4*(mu2*p1 + mu1*(mu2 + p2))))*tau)*(-1 + (-mu1 + mu2 + p1 + p2)/
                                                                                                      sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                         (mu1 - mu2 + p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                         (2 - (mu2*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                 4*(mu2*p1 + mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                            (1/(N*p1))*((-mu1)*n1m1 + mu2*n1m1 - n1m1*p1 - 2*n2m1*p1 + n1m1*p2 + 
                                          n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))*tau)/
                      sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))/
               (2*(4 - (4*mu2*p1)/(mu2*p1 + mu1*(mu2 + p2)) - (4*mu1*p2)/(mu2*p1 + mu1*(mu2 + p2)) - 
                     (2*p1*(-2 + (mu2*(mu1 + mu2 + p1 + p2 - sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                    4*(mu2*p1 + mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                              (1/(N*p1))*(mu1*n1m1 - mu2*n1m1 + n1m1*p1 + 2*n2m1*p1 - n1m1*p2 + 
                                            n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))))/
                     exp((1/2)*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                             4*(mu2*p1 + mu1*(mu2 + p2))))*tau)/sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                                                       4*(mu2*p1 + mu1*(mu2 + p2))) + 
                     ((mu1 - mu2 + p1 - p2 - sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                        (-2 + (mu2*(mu1 + mu2 + p1 + p2 - sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + 
                                                                                              mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                           (1/(N*p1))*(mu1*n1m1 - mu2*n1m1 + n1m1*p1 + 2*n2m1*p1 - n1m1*p2 + 
                                         n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))))/
                     exp((1/2)*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                             4*(mu2*p1 + mu1*(mu2 + p2))))*tau)/sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                                                       4*(mu2*p1 + mu1*(mu2 + p2))) - 
                     (2*exp((1/2)*(-mu1 - mu2 - p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                 4*(mu2*p1 + mu1*(mu2 + p2))))*tau)*p1*
                        (2 - (mu2*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                4*(mu2*p1 + mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                           (1/(N*p1))*((-mu1)*n1m1 + mu2*n1m1 - n1m1*p1 - 2*n2m1*p1 + n1m1*p2 + 
                                         n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))))/
                     sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                     (exp((1/2)*(-mu1 - mu2 - p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                               4*(mu2*p1 + mu1*(mu2 + p2))))*tau)*(mu1 - mu2 + p1 - p2 + 
                                                                                                     sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                        (2 - (mu2*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                4*(mu2*p1 + mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                           (1/(N*p1))*((-mu1)*n1m1 + mu2*n1m1 - n1m1*p1 - 2*n2m1*p1 + n1m1*p2 + 
                                         n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))))/
                     sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))))
  return(gr.p2f)
}

# Gradient equation for mu1
gr.mu1 <- function(par, y){
  p1 <- par[1]
  p2 <- par[2]
  mu1 <- par[3]
  mu2 <- par[4]
  t <- y$t
  tm1 <- c(0, t[1:length(t)-1])
  tau <- t - tm1
  n1 <- y$n1
  n2 <- y$n2
  n3 <- y$n3
  n1m1 <- c(0, n1[1:length(n1)-1])
  n2m1 <- c(0, n2[1:length(n2)-1])
  N <- y$N
  gr.mu1f <- ((n1*(-((2*mu2*(mu2 + p2))/(mu2*p1 + mu1*(mu2 + p2))^2) + 
                     (exp((1/2)*(-mu1 - mu2 - p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                               4*(mu2*p1 + mu1*(mu2 + p2))))*tau)*(mu1 - mu2 + p1 - p2)*
                        (mu1^2*mu2*n1m1 - mu1*mu2^2*n1m1 - mu1*mu2*N*p1 + mu2^2*N*p1 + 2*mu1*mu2*n1m1*p1 - 
                           mu2^2*n1m1*p1 + 2*mu1*mu2*n2m1*p1 - mu2*N*p1^2 + mu2*n1m1*p1^2 + 2*mu2*n2m1*p1^2 + 
                           mu1^2*n1m1*p2 - 2*mu1*mu2*n1m1*p2 - 2*mu1*N*p1*p2 + mu2*N*p1*p2 + mu1*n1m1*p1*p2 - 
                           mu2*n1m1*p1*p2 + 2*mu1*n2m1*p1*p2 - mu1*n1m1*p2^2 - 
                           mu1*mu2*n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                           mu2*N*p1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) - 
                           mu2*n1m1*p1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) - 
                           mu1*n1m1*p2*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))/
                     (N*p1*(mu2*p1 + mu1*(mu2 + p2))*((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))^
                        (3/2)) - ((mu1 - mu2 + p1 - p2)*(mu1^2*mu2*n1m1 - mu1*mu2^2*n1m1 - mu1*mu2*N*p1 + 
                                                           mu2^2*N*p1 + 2*mu1*mu2*n1m1*p1 - mu2^2*n1m1*p1 + 2*mu1*mu2*n2m1*p1 - mu2*N*p1^2 + 
                                                           mu2*n1m1*p1^2 + 2*mu2*n2m1*p1^2 + mu1^2*n1m1*p2 - 2*mu1*mu2*n1m1*p2 - 2*mu1*N*p1*p2 + 
                                                           mu2*N*p1*p2 + mu1*n1m1*p1*p2 - mu2*n1m1*p1*p2 + 2*mu1*n2m1*p1*p2 - mu1*n1m1*p2^2 + 
                                                           mu1*mu2*n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) - 
                                                           mu2*N*p1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                                                           mu2*n1m1*p1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                                                           mu1*n1m1*p2*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))/
                     exp((1/2)*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                             4*(mu2*p1 + mu1*(mu2 + p2))))*tau)/(N*p1*(mu2*p1 + mu1*(mu2 + p2))*
                                                                                                   ((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))^(3/2)) + 
                     (-((mu2*(-1 + (mu1 - mu2 + p1 - p2)/sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + 
                                                                                             mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2))) + 
                        (n1m1*(1 + (mu1 - mu2 + p1 - p2)/sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + 
                                                                                             mu1*(mu2 + p2)))))/(N*p1) - (mu2*(mu2 + p2)*(mu1 + mu2 + p1 + p2 - 
                                                                                                                                            sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))/
                        (mu2*p1 + mu1*(mu2 + p2))^2)/exp((1/2)*(mu1 + mu2 + p1 + p2 + 
                                                                  sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*tau)/
                     sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                     (exp((1/2)*(-mu1 - mu2 - p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                               4*(mu2*p1 + mu1*(mu2 + p2))))*tau)*
                        ((mu2*(-1 + (-mu1 + mu2 - p1 + p2)/sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + 
                                                                                               mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) - 
                           (n1m1*(1 + (-mu1 + mu2 - p1 + p2)/sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + 
                                                                                                 mu1*(mu2 + p2)))))/(N*p1) + (mu2*(mu2 + p2)*(mu1 + mu2 + p1 + p2 + 
                                                                                                                                                sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))/
                           (mu2*p1 + mu1*(mu2 + p2))^2))/sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                4*(mu2*p1 + mu1*(mu2 + p2))) + 
                     ((-1 + (-mu1 + mu2 - p1 + p2)/sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                        (-2 + (mu2*(mu1 + mu2 + p1 + p2 - sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + 
                                                                                              mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                           (1/(N*p1))*(mu1*n1m1 - mu2*n1m1 + n1m1*p1 + 2*n2m1*p1 - n1m1*p2 + 
                                         n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))*tau)/
                     exp((1/2)*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                             4*(mu2*p1 + mu1*(mu2 + p2))))*tau)/
                     (2*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))) + 
                     (exp((1/2)*(-mu1 - mu2 - p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                               4*(mu2*p1 + mu1*(mu2 + p2))))*tau)*(-1 + (mu1 - mu2 + p1 - p2)/
                                                                                                     sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                        (2 - (mu2*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                4*(mu2*p1 + mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                           (1/(N*p1))*((-mu1)*n1m1 + mu2*n1m1 - n1m1*p1 - 2*n2m1*p1 + n1m1*p2 + 
                                         n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))*tau)/
                     (2*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))))/
                ((2*mu2)/(mu2*p1 + mu1*(mu2 + p2)) + 
                   (-2 + (mu2*(mu1 + mu2 + p1 + p2 - sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                            4*(mu2*p1 + mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                      (1/(N*p1))*(mu1*n1m1 - mu2*n1m1 + n1m1*p1 + 2*n2m1*p1 - n1m1*p2 + 
                                    n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))/
                   exp((1/2)*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                           4*(mu2*p1 + mu1*(mu2 + p2))))*tau)/sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                                                     4*(mu2*p1 + mu1*(mu2 + p2))) + 
                   (exp((1/2)*(-mu1 - mu2 - p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                             4*(mu2*p1 + mu1*(mu2 + p2))))*tau)*
                      (2 - (mu2*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                              4*(mu2*p1 + mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                         (1/(N*p1))*((-mu1)*n1m1 + mu2*n1m1 - n1m1*p1 - 2*n2m1*p1 + n1m1*p2 + 
                                       n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))))/
                   sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))) + 
                (n2*(-((8*mu1*p2*(mu2 + p2))/(mu2*p1 + mu1*(mu2 + p2))^2) + (8*p2)/(mu2*p1 + mu1*(mu2 + p2)) - 
                       (2*(mu1 - mu2 + p1 - p2)*(-mu1 + mu2 - p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                               4*(mu2*p1 + mu1*(mu2 + p2))))*(mu1^2*mu2*n1m1 - mu1*mu2^2*n1m1 - mu1*mu2*N*p1 + 
                                                                                                                mu2^2*N*p1 + 2*mu1*mu2*n1m1*p1 - mu2^2*n1m1*p1 + 2*mu1*mu2*n2m1*p1 - mu2*N*p1^2 + 
                                                                                                                mu2*n1m1*p1^2 + 2*mu2*n2m1*p1^2 + mu1^2*n1m1*p2 - 2*mu1*mu2*n1m1*p2 - 2*mu1*N*p1*p2 + 
                                                                                                                mu2*N*p1*p2 + mu1*n1m1*p1*p2 - mu2*n1m1*p1*p2 + 2*mu1*n2m1*p1*p2 - mu1*n1m1*p2^2 + 
                                                                                                                mu1*mu2*n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) - 
                                                                                                                mu2*N*p1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                                                                                                                mu2*n1m1*p1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                                                                                                                mu1*n1m1*p2*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))/
                       exp((1/2)*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                               4*(mu2*p1 + mu1*(mu2 + p2))))*tau)/(N*p1*(mu2*p1 + mu1*(mu2 + p2))*
                                                                                                     ((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))^(3/2)) + 
                       (2*exp((1/2)*(-mu1 - mu2 - p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                   4*(mu2*p1 + mu1*(mu2 + p2))))*tau)*(mu1 - mu2 + p1 - p2)*
                          (mu1 - mu2 + p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                          ((-mu1^2)*mu2*n1m1 + mu1*mu2^2*n1m1 + mu1*mu2*N*p1 - mu2^2*N*p1 - 2*mu1*mu2*n1m1*p1 + 
                             mu2^2*n1m1*p1 - 2*mu1*mu2*n2m1*p1 + mu2*N*p1^2 - mu2*n1m1*p1^2 - 2*mu2*n2m1*p1^2 - 
                             mu1^2*n1m1*p2 + 2*mu1*mu2*n1m1*p2 + 2*mu1*N*p1*p2 - mu2*N*p1*p2 - mu1*n1m1*p1*p2 + 
                             mu2*n1m1*p1*p2 - 2*mu1*n2m1*p1*p2 + mu1*n1m1*p2^2 + 
                             mu1*mu2*n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) - 
                             mu2*N*p1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                             mu2*n1m1*p1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                             mu1*n1m1*p2*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))/
                       (N*p1*(mu2*p1 + mu1*(mu2 + p2))*((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))^
                          (3/2)) - (2*(mu1 - mu2 + p1 - p2 - sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                    4*(mu2*p1 + mu1*(mu2 + p2))))*
                                      (-((mu2*(-1 + (mu1 - mu2 + p1 - p2)/sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                                 4*(mu2*p1 + mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2))) + 
                                         (n1m1*(1 + (mu1 - mu2 + p1 - p2)/sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + 
                                                                                                              mu1*(mu2 + p2)))))/(N*p1) - (mu2*(mu2 + p2)*(mu1 + mu2 + p1 + p2 - 
                                                                                                                                                             sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))/
                                         (mu2*p1 + mu1*(mu2 + p2))^2))/exp((1/2)*(mu1 + mu2 + p1 + p2 + 
                                                                                    sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*tau)/
                       sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) - 
                       (2*exp((1/2)*(-mu1 - mu2 - p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                   4*(mu2*p1 + mu1*(mu2 + p2))))*tau)*(mu1 - mu2 + p1 - p2 + 
                                                                                                         sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                          ((mu2*(-1 + (-mu1 + mu2 - p1 + p2)/sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + 
                                                                                                 mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) - 
                             (n1m1*(1 + (-mu1 + mu2 - p1 + p2)/sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + 
                                                                                                   mu1*(mu2 + p2)))))/(N*p1) + (mu2*(mu2 + p2)*(mu1 + mu2 + p1 + p2 + 
                                                                                                                                                  sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))/
                             (mu2*p1 + mu1*(mu2 + p2))^2))/sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                  4*(mu2*p1 + mu1*(mu2 + p2))) - 
                       (2*(1 + (-mu1 + mu2 - p1 + p2)/sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                          (-2 + (mu2*(mu1 + mu2 + p1 + p2 - sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + 
                                                                                                mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                             (1/(N*p1))*(mu1*n1m1 - mu2*n1m1 + n1m1*p1 + 2*n2m1*p1 - n1m1*p2 + 
                                           n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))))/
                       exp((1/2)*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                               4*(mu2*p1 + mu1*(mu2 + p2))))*tau)/sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                                                         4*(mu2*p1 + mu1*(mu2 + p2))) - 
                       (2*exp((1/2)*(-mu1 - mu2 - p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                   4*(mu2*p1 + mu1*(mu2 + p2))))*tau)*(1 + (mu1 - mu2 + p1 - p2)/
                                                                                                         sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                          (2 - (mu2*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                  4*(mu2*p1 + mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                             (1/(N*p1))*((-mu1)*n1m1 + mu2*n1m1 - n1m1*p1 - 2*n2m1*p1 + n1m1*p2 + 
                                           n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))))/
                       sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) - 
                       ((-1 + (-mu1 + mu2 - p1 + p2)/sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                          (mu1 - mu2 + p1 - p2 - sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                          (-2 + (mu2*(mu1 + mu2 + p1 + p2 - sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + 
                                                                                                mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                             (1/(N*p1))*(mu1*n1m1 - mu2*n1m1 + n1m1*p1 + 2*n2m1*p1 - n1m1*p2 + 
                                           n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))*tau)/
                       exp((1/2)*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                               4*(mu2*p1 + mu1*(mu2 + p2))))*tau)/sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                                                         4*(mu2*p1 + mu1*(mu2 + p2))) - 
                       (exp((1/2)*(-mu1 - mu2 - p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                 4*(mu2*p1 + mu1*(mu2 + p2))))*tau)*(-1 + (mu1 - mu2 + p1 - p2)/
                                                                                                       sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                          (mu1 - mu2 + p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                          (2 - (mu2*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                  4*(mu2*p1 + mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                             (1/(N*p1))*((-mu1)*n1m1 + mu2*n1m1 - n1m1*p1 - 2*n2m1*p1 + n1m1*p2 + 
                                           n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))*tau)/
                       sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))/
                (8*((mu1*p2)/(mu2*p1 + mu1*(mu2 + p2)) - 
                      ((mu1 - mu2 + p1 - p2 - sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                         (-2 + (mu2*(mu1 + mu2 + p1 + p2 - sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + 
                                                                                               mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                            (1/(N*p1))*(mu1*n1m1 - mu2*n1m1 + n1m1*p1 + 2*n2m1*p1 - n1m1*p2 + 
                                          n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))))/
                      exp((1/2)*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                              4*(mu2*p1 + mu1*(mu2 + p2))))*tau)/
                      (4*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))) - 
                      (exp((1/2)*(-mu1 - mu2 - p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                4*(mu2*p1 + mu1*(mu2 + p2))))*tau)*(mu1 - mu2 + p1 - p2 + 
                                                                                                      sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                         (2 - (mu2*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                 4*(mu2*p1 + mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                            (1/(N*p1))*((-mu1)*n1m1 + mu2*n1m1 - n1m1*p1 - 2*n2m1*p1 + n1m1*p2 + 
                                          n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))))/
                      (4*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))) + 
                (n3*((8*mu2*p1*(mu2 + p2))/(mu2*p1 + mu1*(mu2 + p2))^2 + (8*mu1*p2*(mu2 + p2))/
                       (mu2*p1 + mu1*(mu2 + p2))^2 - (8*p2)/(mu2*p1 + mu1*(mu2 + p2)) + 
                       (4*(mu1 - mu2 + p1 - p2)*(mu1^2*mu2*n1m1 - mu1*mu2^2*n1m1 - mu1*mu2*N*p1 + mu2^2*N*p1 + 
                                                   2*mu1*mu2*n1m1*p1 - mu2^2*n1m1*p1 + 2*mu1*mu2*n2m1*p1 - mu2*N*p1^2 + mu2*n1m1*p1^2 + 
                                                   2*mu2*n2m1*p1^2 + mu1^2*n1m1*p2 - 2*mu1*mu2*n1m1*p2 - 2*mu1*N*p1*p2 + mu2*N*p1*p2 + 
                                                   mu1*n1m1*p1*p2 - mu2*n1m1*p1*p2 + 2*mu1*n2m1*p1*p2 - mu1*n1m1*p2^2 + 
                                                   mu1*mu2*n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) - 
                                                   mu2*N*p1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                                                   mu2*n1m1*p1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                                                   mu1*n1m1*p2*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))/
                       exp((1/2)*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                               4*(mu2*p1 + mu1*(mu2 + p2))))*tau)/(N*(mu2*p1 + mu1*(mu2 + p2))*
                                                                                                     ((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))^(3/2)) + 
                       (2*(mu1 - mu2 + p1 - p2)*(-mu1 + mu2 - p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                               4*(mu2*p1 + mu1*(mu2 + p2))))*(mu1^2*mu2*n1m1 - mu1*mu2^2*n1m1 - mu1*mu2*N*p1 + 
                                                                                                                mu2^2*N*p1 + 2*mu1*mu2*n1m1*p1 - mu2^2*n1m1*p1 + 2*mu1*mu2*n2m1*p1 - mu2*N*p1^2 + 
                                                                                                                mu2*n1m1*p1^2 + 2*mu2*n2m1*p1^2 + mu1^2*n1m1*p2 - 2*mu1*mu2*n1m1*p2 - 2*mu1*N*p1*p2 + 
                                                                                                                mu2*N*p1*p2 + mu1*n1m1*p1*p2 - mu2*n1m1*p1*p2 + 2*mu1*n2m1*p1*p2 - mu1*n1m1*p2^2 + 
                                                                                                                mu1*mu2*n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) - 
                                                                                                                mu2*N*p1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                                                                                                                mu2*n1m1*p1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                                                                                                                mu1*n1m1*p2*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))/
                       exp((1/2)*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                               4*(mu2*p1 + mu1*(mu2 + p2))))*tau)/(N*p1*(mu2*p1 + mu1*(mu2 + p2))*
                                                                                                     ((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))^(3/2)) + 
                       (4*exp((1/2)*(-mu1 - mu2 - p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                   4*(mu2*p1 + mu1*(mu2 + p2))))*tau)*(mu1 - mu2 + p1 - p2)*
                          ((-mu1^2)*mu2*n1m1 + mu1*mu2^2*n1m1 + mu1*mu2*N*p1 - mu2^2*N*p1 - 2*mu1*mu2*n1m1*p1 + 
                             mu2^2*n1m1*p1 - 2*mu1*mu2*n2m1*p1 + mu2*N*p1^2 - mu2*n1m1*p1^2 - 2*mu2*n2m1*p1^2 - 
                             mu1^2*n1m1*p2 + 2*mu1*mu2*n1m1*p2 + 2*mu1*N*p1*p2 - mu2*N*p1*p2 - mu1*n1m1*p1*p2 + 
                             mu2*n1m1*p1*p2 - 2*mu1*n2m1*p1*p2 + mu1*n1m1*p2^2 + 
                             mu1*mu2*n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) - 
                             mu2*N*p1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                             mu2*n1m1*p1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                             mu1*n1m1*p2*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))/
                       (N*(mu2*p1 + mu1*(mu2 + p2))*((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))^
                          (3/2)) - (2*exp((1/2)*(-mu1 - mu2 - p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                               4*(mu2*p1 + mu1*(mu2 + p2))))*tau)*(mu1 - mu2 + p1 - p2)*
                                      (mu1 - mu2 + p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                                      ((-mu1^2)*mu2*n1m1 + mu1*mu2^2*n1m1 + mu1*mu2*N*p1 - mu2^2*N*p1 - 2*mu1*mu2*n1m1*p1 + 
                                         mu2^2*n1m1*p1 - 2*mu1*mu2*n2m1*p1 + mu2*N*p1^2 - mu2*n1m1*p1^2 - 2*mu2*n2m1*p1^2 - 
                                         mu1^2*n1m1*p2 + 2*mu1*mu2*n1m1*p2 + 2*mu1*N*p1*p2 - mu2*N*p1*p2 - mu1*n1m1*p1*p2 + 
                                         mu2*n1m1*p1*p2 - 2*mu1*n2m1*p1*p2 + mu1*n1m1*p2^2 + 
                                         mu1*mu2*n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) - 
                                         mu2*N*p1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                                         mu2*n1m1*p1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                                         mu1*n1m1*p2*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))/
                       (N*p1*(mu2*p1 + mu1*(mu2 + p2))*((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))^
                          (3/2)) - 
                       (4*p1*(-((mu2*(-1 + (mu1 - mu2 + p1 - p2)/sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                        4*(mu2*p1 + mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2))) + 
                                (n1m1*(1 + (mu1 - mu2 + p1 - p2)/sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + 
                                                                                                     mu1*(mu2 + p2)))))/(N*p1) - (mu2*(mu2 + p2)*(mu1 + mu2 + p1 + p2 - 
                                                                                                                                                    sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))/
                                (mu2*p1 + mu1*(mu2 + p2))^2))/exp((1/2)*(mu1 + mu2 + p1 + p2 + 
                                                                           sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*tau)/
                       sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                       (2*(mu1 - mu2 + p1 - p2 - sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                          (-((mu2*(-1 + (mu1 - mu2 + p1 - p2)/sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                     4*(mu2*p1 + mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2))) + 
                             (n1m1*(1 + (mu1 - mu2 + p1 - p2)/sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + 
                                                                                                  mu1*(mu2 + p2)))))/(N*p1) - (mu2*(mu2 + p2)*(mu1 + mu2 + p1 + p2 - 
                                                                                                                                                 sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))/
                             (mu2*p1 + mu1*(mu2 + p2))^2))/exp((1/2)*(mu1 + mu2 + p1 + p2 + 
                                                                        sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*tau)/
                       sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) - 
                       (4*exp((1/2)*(-mu1 - mu2 - p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                   4*(mu2*p1 + mu1*(mu2 + p2))))*tau)*p1*
                          ((mu2*(-1 + (-mu1 + mu2 - p1 + p2)/sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + 
                                                                                                 mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) - 
                             (n1m1*(1 + (-mu1 + mu2 - p1 + p2)/sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + 
                                                                                                   mu1*(mu2 + p2)))))/(N*p1) + (mu2*(mu2 + p2)*(mu1 + mu2 + p1 + p2 + 
                                                                                                                                                  sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))/
                             (mu2*p1 + mu1*(mu2 + p2))^2))/sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                  4*(mu2*p1 + mu1*(mu2 + p2))) + 
                       (2*exp((1/2)*(-mu1 - mu2 - p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                   4*(mu2*p1 + mu1*(mu2 + p2))))*tau)*(mu1 - mu2 + p1 - p2 + 
                                                                                                         sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                          ((mu2*(-1 + (-mu1 + mu2 - p1 + p2)/sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + 
                                                                                                 mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) - 
                             (n1m1*(1 + (-mu1 + mu2 - p1 + p2)/sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + 
                                                                                                   mu1*(mu2 + p2)))))/(N*p1) + (mu2*(mu2 + p2)*(mu1 + mu2 + p1 + p2 + 
                                                                                                                                                  sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))/
                             (mu2*p1 + mu1*(mu2 + p2))^2))/sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                  4*(mu2*p1 + mu1*(mu2 + p2))) + 
                       (2*(1 + (-mu1 + mu2 - p1 + p2)/sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                          (-2 + (mu2*(mu1 + mu2 + p1 + p2 - sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + 
                                                                                                mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                             (1/(N*p1))*(mu1*n1m1 - mu2*n1m1 + n1m1*p1 + 2*n2m1*p1 - n1m1*p2 + 
                                           n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))))/
                       exp((1/2)*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                               4*(mu2*p1 + mu1*(mu2 + p2))))*tau)/sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                                                         4*(mu2*p1 + mu1*(mu2 + p2))) + 
                       (2*exp((1/2)*(-mu1 - mu2 - p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                   4*(mu2*p1 + mu1*(mu2 + p2))))*tau)*(1 + (mu1 - mu2 + p1 - p2)/
                                                                                                         sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                          (2 - (mu2*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                  4*(mu2*p1 + mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                             (1/(N*p1))*((-mu1)*n1m1 + mu2*n1m1 - n1m1*p1 - 2*n2m1*p1 + n1m1*p2 + 
                                           n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))))/
                       sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) - 
                       (2*p1*(-1 + (-mu1 + mu2 - p1 + p2)/sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                 4*(mu2*p1 + mu1*(mu2 + p2))))*
                          (-2 + (mu2*(mu1 + mu2 + p1 + p2 - sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + 
                                                                                                mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                             (1/(N*p1))*(mu1*n1m1 - mu2*n1m1 + n1m1*p1 + 2*n2m1*p1 - n1m1*p2 + 
                                           n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))*tau)/
                       exp((1/2)*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                               4*(mu2*p1 + mu1*(mu2 + p2))))*tau)/sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                                                         4*(mu2*p1 + mu1*(mu2 + p2))) + 
                       ((-1 + (-mu1 + mu2 - p1 + p2)/sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                          (mu1 - mu2 + p1 - p2 - sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                          (-2 + (mu2*(mu1 + mu2 + p1 + p2 - sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + 
                                                                                                mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                             (1/(N*p1))*(mu1*n1m1 - mu2*n1m1 + n1m1*p1 + 2*n2m1*p1 - n1m1*p2 + 
                                           n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))*tau)/
                       exp((1/2)*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                               4*(mu2*p1 + mu1*(mu2 + p2))))*tau)/sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                                                         4*(mu2*p1 + mu1*(mu2 + p2))) - 
                       (2*exp((1/2)*(-mu1 - mu2 - p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                   4*(mu2*p1 + mu1*(mu2 + p2))))*tau)*p1*(-1 + (mu1 - mu2 + p1 - p2)/
                                                                                                            sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                          (2 - (mu2*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                  4*(mu2*p1 + mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                             (1/(N*p1))*((-mu1)*n1m1 + mu2*n1m1 - n1m1*p1 - 2*n2m1*p1 + n1m1*p2 + 
                                           n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))*tau)/
                       sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                       (exp((1/2)*(-mu1 - mu2 - p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                 4*(mu2*p1 + mu1*(mu2 + p2))))*tau)*(-1 + (mu1 - mu2 + p1 - p2)/
                                                                                                       sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                          (mu1 - mu2 + p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                          (2 - (mu2*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                  4*(mu2*p1 + mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                             (1/(N*p1))*((-mu1)*n1m1 + mu2*n1m1 - n1m1*p1 - 2*n2m1*p1 + n1m1*p2 + 
                                           n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))*tau)/
                       sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))/
                (2*(4 - (4*mu2*p1)/(mu2*p1 + mu1*(mu2 + p2)) - (4*mu1*p2)/(mu2*p1 + mu1*(mu2 + p2)) - 
                      (2*p1*(-2 + (mu2*(mu1 + mu2 + p1 + p2 - sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                     4*(mu2*p1 + mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                               (1/(N*p1))*(mu1*n1m1 - mu2*n1m1 + n1m1*p1 + 2*n2m1*p1 - n1m1*p2 + 
                                             n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))))/
                      exp((1/2)*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                              4*(mu2*p1 + mu1*(mu2 + p2))))*tau)/sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                                                        4*(mu2*p1 + mu1*(mu2 + p2))) + 
                      ((mu1 - mu2 + p1 - p2 - sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                         (-2 + (mu2*(mu1 + mu2 + p1 + p2 - sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + 
                                                                                               mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                            (1/(N*p1))*(mu1*n1m1 - mu2*n1m1 + n1m1*p1 + 2*n2m1*p1 - n1m1*p2 + 
                                          n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))))/
                      exp((1/2)*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                              4*(mu2*p1 + mu1*(mu2 + p2))))*tau)/sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                                                        4*(mu2*p1 + mu1*(mu2 + p2))) - 
                      (2*exp((1/2)*(-mu1 - mu2 - p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                  4*(mu2*p1 + mu1*(mu2 + p2))))*tau)*p1*
                         (2 - (mu2*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                 4*(mu2*p1 + mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                            (1/(N*p1))*((-mu1)*n1m1 + mu2*n1m1 - n1m1*p1 - 2*n2m1*p1 + n1m1*p2 + 
                                          n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))))/
                      sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                      (exp((1/2)*(-mu1 - mu2 - p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                4*(mu2*p1 + mu1*(mu2 + p2))))*tau)*(mu1 - mu2 + p1 - p2 + 
                                                                                                      sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                         (2 - (mu2*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                 4*(mu2*p1 + mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                            (1/(N*p1))*((-mu1)*n1m1 + mu2*n1m1 - n1m1*p1 - 2*n2m1*p1 + n1m1*p2 + 
                                          n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))))/
                      sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))))
  return(gr.mu1f)
}

# Gradient equation for mu2
gr.mu2 <- function(par, y){
  p1 <- par[1]
  p2 <- par[2]
  mu1 <- par[3]
  mu2 <- par[4]
  t <- y$t
  tm1 <- c(0, t[1:length(t)-1])
  tau <- t - tm1
  n1 <- y$n1
  n2 <- y$n2
  n3 <- y$n3
  n1m1 <- c(0, n1[1:length(n1)-1])
  n2m1 <- c(0, n2[1:length(n2)-1])
  N <- y$N
  gr.mu2f <- ((n2*(-((8*mu1*(mu1 + p1)*p2)/(mu2*p1 + mu1*(mu2 + p2))^2) + 
                     (2*(mu1 - mu2 + p1 - p2)*(-mu1 + mu2 - p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                             4*(mu2*p1 + mu1*(mu2 + p2))))*(mu1^2*mu2*n1m1 - mu1*mu2^2*n1m1 - mu1*mu2*N*p1 + 
                                                                                                              mu2^2*N*p1 + 2*mu1*mu2*n1m1*p1 - mu2^2*n1m1*p1 + 2*mu1*mu2*n2m1*p1 - mu2*N*p1^2 + 
                                                                                                              mu2*n1m1*p1^2 + 2*mu2*n2m1*p1^2 + mu1^2*n1m1*p2 - 2*mu1*mu2*n1m1*p2 - 2*mu1*N*p1*p2 + 
                                                                                                              mu2*N*p1*p2 + mu1*n1m1*p1*p2 - mu2*n1m1*p1*p2 + 2*mu1*n2m1*p1*p2 - mu1*n1m1*p2^2 + 
                                                                                                              mu1*mu2*n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) - 
                                                                                                              mu2*N*p1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                                                                                                              mu2*n1m1*p1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                                                                                                              mu1*n1m1*p2*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))/
                     exp((1/2)*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                             4*(mu2*p1 + mu1*(mu2 + p2))))*tau)/(N*p1*(mu2*p1 + mu1*(mu2 + p2))*
                                                                                                   ((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))^(3/2)) - 
                     (2*exp((1/2)*(-mu1 - mu2 - p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                 4*(mu2*p1 + mu1*(mu2 + p2))))*tau)*(mu1 - mu2 + p1 - p2)*
                        (mu1 - mu2 + p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                        ((-mu1^2)*mu2*n1m1 + mu1*mu2^2*n1m1 + mu1*mu2*N*p1 - mu2^2*N*p1 - 2*mu1*mu2*n1m1*p1 + 
                           mu2^2*n1m1*p1 - 2*mu1*mu2*n2m1*p1 + mu2*N*p1^2 - mu2*n1m1*p1^2 - 2*mu2*n2m1*p1^2 - 
                           mu1^2*n1m1*p2 + 2*mu1*mu2*n1m1*p2 + 2*mu1*N*p1*p2 - mu2*N*p1*p2 - mu1*n1m1*p1*p2 + 
                           mu2*n1m1*p1*p2 - 2*mu1*n2m1*p1*p2 + mu1*n1m1*p2^2 + 
                           mu1*mu2*n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) - 
                           mu2*N*p1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                           mu2*n1m1*p1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                           mu1*n1m1*p2*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))/
                     (N*p1*(mu2*p1 + mu1*(mu2 + p2))*((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))^
                        (3/2)) - (2*(mu1 - mu2 + p1 - p2 - sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                  4*(mu2*p1 + mu1*(mu2 + p2))))*
                                    ((n1m1*(-1 + (-mu1 + mu2 - p1 + p2)/sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + 
                                                                                                            mu1*(mu2 + p2)))))/(N*p1) - 
                                       (mu2*(-1 + (-mu1 + mu2 - p1 + p2)/sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + 
                                                                                                             mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) - 
                                       (mu2*(mu1 + p1)*(mu1 + mu2 + p1 + p2 - sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                                     4*(mu2*p1 + mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2))^2 + 
                                       (mu1 + mu2 + p1 + p2 - sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))/
                                       (mu2*p1 + mu1*(mu2 + p2))))/exp((1/2)*(mu1 + mu2 + p1 + p2 + 
                                                                                sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*tau)/
                     sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) - 
                     (2*exp((1/2)*(-mu1 - mu2 - p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                 4*(mu2*p1 + mu1*(mu2 + p2))))*tau)*(mu1 - mu2 + p1 - p2 + 
                                                                                                       sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                        (-((n1m1*(-1 + (mu1 - mu2 + p1 - p2)/sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + 
                                                                                                 mu1*(mu2 + p2)))))/(N*p1)) + 
                           (mu2*(-1 + (mu1 - mu2 + p1 - p2)/sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + 
                                                                                                mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                           (mu2*(mu1 + p1)*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                         4*(mu2*p1 + mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2))^2 - 
                           (mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))/
                           (mu2*p1 + mu1*(mu2 + p2))))/sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                              4*(mu2*p1 + mu1*(mu2 + p2))) - 
                     (2*(-1 + (mu1 - mu2 + p1 - p2)/sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                        (-2 + (mu2*(mu1 + mu2 + p1 + p2 - sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + 
                                                                                              mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                           (1/(N*p1))*(mu1*n1m1 - mu2*n1m1 + n1m1*p1 + 2*n2m1*p1 - n1m1*p2 + 
                                         n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))))/
                     exp((1/2)*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                             4*(mu2*p1 + mu1*(mu2 + p2))))*tau)/sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                                                       4*(mu2*p1 + mu1*(mu2 + p2))) - 
                     (2*exp((1/2)*(-mu1 - mu2 - p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                 4*(mu2*p1 + mu1*(mu2 + p2))))*tau)*(-1 + (-mu1 + mu2 - p1 + p2)/
                                                                                                       sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                        (2 - (mu2*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                4*(mu2*p1 + mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                           (1/(N*p1))*((-mu1)*n1m1 + mu2*n1m1 - n1m1*p1 - 2*n2m1*p1 + n1m1*p2 + 
                                         n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))))/
                     sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) - 
                     ((-mu1 + mu2 - p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))^2*
                        (mu1^2*mu2*n1m1 - mu1*mu2^2*n1m1 - mu1*mu2*N*p1 + mu2^2*N*p1 + 2*mu1*mu2*n1m1*p1 - 
                           mu2^2*n1m1*p1 + 2*mu1*mu2*n2m1*p1 - mu2*N*p1^2 + mu2*n1m1*p1^2 + 2*mu2*n2m1*p1^2 + 
                           mu1^2*n1m1*p2 - 2*mu1*mu2*n1m1*p2 - 2*mu1*N*p1*p2 + mu2*N*p1*p2 + mu1*n1m1*p1*p2 - 
                           mu2*n1m1*p1*p2 + 2*mu1*n2m1*p1*p2 - mu1*n1m1*p2^2 + mu1*mu2*n1m1*
                           sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) - 
                           mu2*N*p1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                           mu2*n1m1*p1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                           mu1*n1m1*p2*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*tau)/
                     exp((1/2)*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                             4*(mu2*p1 + mu1*(mu2 + p2))))*tau)/(N*p1*(mu2*p1 + mu1*(mu2 + p2))*
                                                                                                   ((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))) + 
                     (exp((1/2)*(-mu1 - mu2 - p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                               4*(mu2*p1 + mu1*(mu2 + p2))))*tau)*(mu1 - mu2 + p1 - p2 + 
                                                                                                     sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))^2*
                        ((-mu1^2)*mu2*n1m1 + mu1*mu2^2*n1m1 + mu1*mu2*N*p1 - mu2^2*N*p1 - 2*mu1*mu2*n1m1*p1 + 
                           mu2^2*n1m1*p1 - 2*mu1*mu2*n2m1*p1 + mu2*N*p1^2 - mu2*n1m1*p1^2 - 2*mu2*n2m1*p1^2 - 
                           mu1^2*n1m1*p2 + 2*mu1*mu2*n1m1*p2 + 2*mu1*N*p1*p2 - mu2*N*p1*p2 - mu1*n1m1*p1*p2 + 
                           mu2*n1m1*p1*p2 - 2*mu1*n2m1*p1*p2 + mu1*n1m1*p2^2 + 
                           mu1*mu2*n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) - 
                           mu2*N*p1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                           mu2*n1m1*p1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                           mu1*n1m1*p2*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*tau)/
                     (N*p1*(mu2*p1 + mu1*(mu2 + p2))*((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))))/
                (8*((mu1*p2)/(mu2*p1 + mu1*(mu2 + p2)) - 
                      ((mu1 - mu2 + p1 - p2 - sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                         (-2 + (mu2*(mu1 + mu2 + p1 + p2 - sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + 
                                                                                               mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                            (1/(N*p1))*(mu1*n1m1 - mu2*n1m1 + n1m1*p1 + 2*n2m1*p1 - n1m1*p2 + 
                                          n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))))/
                      exp((1/2)*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                              4*(mu2*p1 + mu1*(mu2 + p2))))*tau)/
                      (4*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))) - 
                      (exp((1/2)*(-mu1 - mu2 - p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                4*(mu2*p1 + mu1*(mu2 + p2))))*tau)*(mu1 - mu2 + p1 - p2 + 
                                                                                                      sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                         (2 - (mu2*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                 4*(mu2*p1 + mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                            (1/(N*p1))*((-mu1)*n1m1 + mu2*n1m1 - n1m1*p1 - 2*n2m1*p1 + n1m1*p2 + 
                                          n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))))/
                      (4*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))) + 
                (n1*(-((2*mu2*(mu1 + p1))/(mu2*p1 + mu1*(mu2 + p2))^2) + 2/(mu2*p1 + mu1*(mu2 + p2)) + 
                       ((mu1 - mu2 + p1 - p2)*(mu1^2*mu2*n1m1 - mu1*mu2^2*n1m1 - mu1*mu2*N*p1 + mu2^2*N*p1 + 
                                                 2*mu1*mu2*n1m1*p1 - mu2^2*n1m1*p1 + 2*mu1*mu2*n2m1*p1 - mu2*N*p1^2 + mu2*n1m1*p1^2 + 
                                                 2*mu2*n2m1*p1^2 + mu1^2*n1m1*p2 - 2*mu1*mu2*n1m1*p2 - 2*mu1*N*p1*p2 + mu2*N*p1*p2 + 
                                                 mu1*n1m1*p1*p2 - mu2*n1m1*p1*p2 + 2*mu1*n2m1*p1*p2 - mu1*n1m1*p2^2 + 
                                                 mu1*mu2*n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) - 
                                                 mu2*N*p1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                                                 mu2*n1m1*p1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                                                 mu1*n1m1*p2*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))/
                       exp((1/2)*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                               4*(mu2*p1 + mu1*(mu2 + p2))))*tau)/(N*p1*(mu2*p1 + mu1*(mu2 + p2))*
                                                                                                     ((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))^(3/2)) + 
                       (exp((1/2)*(-mu1 - mu2 - p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                 4*(mu2*p1 + mu1*(mu2 + p2))))*tau)*(mu1 - mu2 + p1 - p2)*
                          ((-mu1^2)*mu2*n1m1 + mu1*mu2^2*n1m1 + mu1*mu2*N*p1 - mu2^2*N*p1 - 2*mu1*mu2*n1m1*p1 + 
                             mu2^2*n1m1*p1 - 2*mu1*mu2*n2m1*p1 + mu2*N*p1^2 - mu2*n1m1*p1^2 - 2*mu2*n2m1*p1^2 - 
                             mu1^2*n1m1*p2 + 2*mu1*mu2*n1m1*p2 + 2*mu1*N*p1*p2 - mu2*N*p1*p2 - mu1*n1m1*p1*p2 + 
                             mu2*n1m1*p1*p2 - 2*mu1*n2m1*p1*p2 + mu1*n1m1*p2^2 + 
                             mu1*mu2*n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) - 
                             mu2*N*p1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                             mu2*n1m1*p1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                             mu1*n1m1*p2*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))/
                       (N*p1*(mu2*p1 + mu1*(mu2 + p2))*((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))^
                          (3/2)) + 
                       ((n1m1*(-1 + (-mu1 + mu2 - p1 + p2)/sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + 
                                                                                               mu1*(mu2 + p2)))))/(N*p1) - 
                          (mu2*(-1 + (-mu1 + mu2 - p1 + p2)/sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + 
                                                                                                mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) - 
                          (mu2*(mu1 + p1)*(mu1 + mu2 + p1 + p2 - sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                        4*(mu2*p1 + mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2))^2 + 
                          (mu1 + mu2 + p1 + p2 - sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))/
                          (mu2*p1 + mu1*(mu2 + p2)))/exp((1/2)*(mu1 + mu2 + p1 + p2 + 
                                                                  sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*tau)/
                       sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                       (exp((1/2)*(-mu1 - mu2 - p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                 4*(mu2*p1 + mu1*(mu2 + p2))))*tau)*
                          (-((n1m1*(-1 + (mu1 - mu2 + p1 - p2)/sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + 
                                                                                                   mu1*(mu2 + p2)))))/(N*p1)) + 
                             (mu2*(-1 + (mu1 - mu2 + p1 - p2)/sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + 
                                                                                                  mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                             (mu2*(mu1 + p1)*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                           4*(mu2*p1 + mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2))^2 - 
                             (mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))/
                             (mu2*p1 + mu1*(mu2 + p2))))/sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                4*(mu2*p1 + mu1*(mu2 + p2))) + 
                       ((-1 + (mu1 - mu2 + p1 - p2)/sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                          (-2 + (mu2*(mu1 + mu2 + p1 + p2 - sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + 
                                                                                                mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                             (1/(N*p1))*(mu1*n1m1 - mu2*n1m1 + n1m1*p1 + 2*n2m1*p1 - n1m1*p2 + 
                                           n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))*tau)/
                       exp((1/2)*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                               4*(mu2*p1 + mu1*(mu2 + p2))))*tau)/
                       (2*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))) + 
                       (exp((1/2)*(-mu1 - mu2 - p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                 4*(mu2*p1 + mu1*(mu2 + p2))))*tau)*(-1 + (-mu1 + mu2 - p1 + p2)/
                                                                                                       sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                          (2 - (mu2*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                  4*(mu2*p1 + mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                             (1/(N*p1))*((-mu1)*n1m1 + mu2*n1m1 - n1m1*p1 - 2*n2m1*p1 + n1m1*p2 + 
                                           n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))*tau)/
                       (2*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))))/
                ((2*mu2)/(mu2*p1 + mu1*(mu2 + p2)) + 
                   (-2 + (mu2*(mu1 + mu2 + p1 + p2 - sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                            4*(mu2*p1 + mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                      (1/(N*p1))*(mu1*n1m1 - mu2*n1m1 + n1m1*p1 + 2*n2m1*p1 - n1m1*p2 + 
                                    n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))/
                   exp((1/2)*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                           4*(mu2*p1 + mu1*(mu2 + p2))))*tau)/sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                                                     4*(mu2*p1 + mu1*(mu2 + p2))) + 
                   (exp((1/2)*(-mu1 - mu2 - p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                             4*(mu2*p1 + mu1*(mu2 + p2))))*tau)*
                      (2 - (mu2*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                              4*(mu2*p1 + mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                         (1/(N*p1))*((-mu1)*n1m1 + mu2*n1m1 - n1m1*p1 - 2*n2m1*p1 + n1m1*p2 + 
                                       n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))))/
                   sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))) + 
                (n3*((8*mu2*p1*(mu1 + p1))/(mu2*p1 + mu1*(mu2 + p2))^2 + (8*mu1*(mu1 + p1)*p2)/
                       (mu2*p1 + mu1*(mu2 + p2))^2 - (8*p1)/(mu2*p1 + mu1*(mu2 + p2)) - 
                       (4*(mu1 - mu2 + p1 - p2)*(mu1^2*mu2*n1m1 - mu1*mu2^2*n1m1 - mu1*mu2*N*p1 + mu2^2*N*p1 + 
                                                   2*mu1*mu2*n1m1*p1 - mu2^2*n1m1*p1 + 2*mu1*mu2*n2m1*p1 - mu2*N*p1^2 + mu2*n1m1*p1^2 + 
                                                   2*mu2*n2m1*p1^2 + mu1^2*n1m1*p2 - 2*mu1*mu2*n1m1*p2 - 2*mu1*N*p1*p2 + mu2*N*p1*p2 + 
                                                   mu1*n1m1*p1*p2 - mu2*n1m1*p1*p2 + 2*mu1*n2m1*p1*p2 - mu1*n1m1*p2^2 + 
                                                   mu1*mu2*n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) - 
                                                   mu2*N*p1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                                                   mu2*n1m1*p1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                                                   mu1*n1m1*p2*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))/
                       exp((1/2)*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                               4*(mu2*p1 + mu1*(mu2 + p2))))*tau)/(N*(mu2*p1 + mu1*(mu2 + p2))*
                                                                                                     ((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))^(3/2)) - 
                       (2*(mu1 - mu2 + p1 - p2)*(-mu1 + mu2 - p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                               4*(mu2*p1 + mu1*(mu2 + p2))))*(mu1^2*mu2*n1m1 - mu1*mu2^2*n1m1 - mu1*mu2*N*p1 + 
                                                                                                                mu2^2*N*p1 + 2*mu1*mu2*n1m1*p1 - mu2^2*n1m1*p1 + 2*mu1*mu2*n2m1*p1 - mu2*N*p1^2 + 
                                                                                                                mu2*n1m1*p1^2 + 2*mu2*n2m1*p1^2 + mu1^2*n1m1*p2 - 2*mu1*mu2*n1m1*p2 - 2*mu1*N*p1*p2 + 
                                                                                                                mu2*N*p1*p2 + mu1*n1m1*p1*p2 - mu2*n1m1*p1*p2 + 2*mu1*n2m1*p1*p2 - mu1*n1m1*p2^2 + 
                                                                                                                mu1*mu2*n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) - 
                                                                                                                mu2*N*p1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                                                                                                                mu2*n1m1*p1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                                                                                                                mu1*n1m1*p2*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))/
                       exp((1/2)*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                               4*(mu2*p1 + mu1*(mu2 + p2))))*tau)/(N*p1*(mu2*p1 + mu1*(mu2 + p2))*
                                                                                                     ((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))^(3/2)) - 
                       (4*exp((1/2)*(-mu1 - mu2 - p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                   4*(mu2*p1 + mu1*(mu2 + p2))))*tau)*(mu1 - mu2 + p1 - p2)*
                          ((-mu1^2)*mu2*n1m1 + mu1*mu2^2*n1m1 + mu1*mu2*N*p1 - mu2^2*N*p1 - 2*mu1*mu2*n1m1*p1 + 
                             mu2^2*n1m1*p1 - 2*mu1*mu2*n2m1*p1 + mu2*N*p1^2 - mu2*n1m1*p1^2 - 2*mu2*n2m1*p1^2 - 
                             mu1^2*n1m1*p2 + 2*mu1*mu2*n1m1*p2 + 2*mu1*N*p1*p2 - mu2*N*p1*p2 - mu1*n1m1*p1*p2 + 
                             mu2*n1m1*p1*p2 - 2*mu1*n2m1*p1*p2 + mu1*n1m1*p2^2 + 
                             mu1*mu2*n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) - 
                             mu2*N*p1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                             mu2*n1m1*p1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                             mu1*n1m1*p2*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))/
                       (N*(mu2*p1 + mu1*(mu2 + p2))*((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))^
                          (3/2)) + (2*exp((1/2)*(-mu1 - mu2 - p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                               4*(mu2*p1 + mu1*(mu2 + p2))))*tau)*(mu1 - mu2 + p1 - p2)*
                                      (mu1 - mu2 + p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                                      ((-mu1^2)*mu2*n1m1 + mu1*mu2^2*n1m1 + mu1*mu2*N*p1 - mu2^2*N*p1 - 2*mu1*mu2*n1m1*p1 + 
                                         mu2^2*n1m1*p1 - 2*mu1*mu2*n2m1*p1 + mu2*N*p1^2 - mu2*n1m1*p1^2 - 2*mu2*n2m1*p1^2 - 
                                         mu1^2*n1m1*p2 + 2*mu1*mu2*n1m1*p2 + 2*mu1*N*p1*p2 - mu2*N*p1*p2 - mu1*n1m1*p1*p2 + 
                                         mu2*n1m1*p1*p2 - 2*mu1*n2m1*p1*p2 + mu1*n1m1*p2^2 + 
                                         mu1*mu2*n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) - 
                                         mu2*N*p1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                                         mu2*n1m1*p1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                                         mu1*n1m1*p2*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))/
                       (N*p1*(mu2*p1 + mu1*(mu2 + p2))*((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))^
                          (3/2)) - 
                       (4*p1*((n1m1*(-1 + (-mu1 + mu2 - p1 + p2)/sqrt((mu1 + mu2 + p1 + p2)^2 - 4*
                                                                        (mu2*p1 + mu1*(mu2 + p2)))))/(N*p1) - 
                                (mu2*(-1 + (-mu1 + mu2 - p1 + p2)/sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + 
                                                                                                      mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) - 
                                (mu2*(mu1 + p1)*(mu1 + mu2 + p1 + p2 - sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                              4*(mu2*p1 + mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2))^2 + 
                                (mu1 + mu2 + p1 + p2 - sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))/
                                (mu2*p1 + mu1*(mu2 + p2))))/exp((1/2)*(mu1 + mu2 + p1 + p2 + 
                                                                         sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*tau)/
                       sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                       (2*(mu1 - mu2 + p1 - p2 - sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                          ((n1m1*(-1 + (-mu1 + mu2 - p1 + p2)/sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + 
                                                                                                  mu1*(mu2 + p2)))))/(N*p1) - 
                             (mu2*(-1 + (-mu1 + mu2 - p1 + p2)/sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + 
                                                                                                   mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) - 
                             (mu2*(mu1 + p1)*(mu1 + mu2 + p1 + p2 - sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                           4*(mu2*p1 + mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2))^2 + 
                             (mu1 + mu2 + p1 + p2 - sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))/
                             (mu2*p1 + mu1*(mu2 + p2))))/exp((1/2)*(mu1 + mu2 + p1 + p2 + 
                                                                      sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*tau)/
                       sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) - 
                       (4*exp((1/2)*(-mu1 - mu2 - p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                   4*(mu2*p1 + mu1*(mu2 + p2))))*tau)*p1*
                          (-((n1m1*(-1 + (mu1 - mu2 + p1 - p2)/sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + 
                                                                                                   mu1*(mu2 + p2)))))/(N*p1)) + 
                             (mu2*(-1 + (mu1 - mu2 + p1 - p2)/sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + 
                                                                                                  mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                             (mu2*(mu1 + p1)*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                           4*(mu2*p1 + mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2))^2 - 
                             (mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))/
                             (mu2*p1 + mu1*(mu2 + p2))))/sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                4*(mu2*p1 + mu1*(mu2 + p2))) + 
                       (2*exp((1/2)*(-mu1 - mu2 - p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                   4*(mu2*p1 + mu1*(mu2 + p2))))*tau)*(mu1 - mu2 + p1 - p2 + 
                                                                                                         sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                          (-((n1m1*(-1 + (mu1 - mu2 + p1 - p2)/sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + 
                                                                                                   mu1*(mu2 + p2)))))/(N*p1)) + 
                             (mu2*(-1 + (mu1 - mu2 + p1 - p2)/sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + 
                                                                                                  mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                             (mu2*(mu1 + p1)*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                           4*(mu2*p1 + mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2))^2 - 
                             (mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))/
                             (mu2*p1 + mu1*(mu2 + p2))))/sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                4*(mu2*p1 + mu1*(mu2 + p2))) + 
                       (2*(-1 + (mu1 - mu2 + p1 - p2)/sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                          (-2 + (mu2*(mu1 + mu2 + p1 + p2 - sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + 
                                                                                                mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                             (1/(N*p1))*(mu1*n1m1 - mu2*n1m1 + n1m1*p1 + 2*n2m1*p1 - n1m1*p2 + 
                                           n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))))/
                       exp((1/2)*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                               4*(mu2*p1 + mu1*(mu2 + p2))))*tau)/sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                                                         4*(mu2*p1 + mu1*(mu2 + p2))) + 
                       (2*exp((1/2)*(-mu1 - mu2 - p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                   4*(mu2*p1 + mu1*(mu2 + p2))))*tau)*(-1 + (-mu1 + mu2 - p1 + p2)/
                                                                                                         sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                          (2 - (mu2*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                  4*(mu2*p1 + mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                             (1/(N*p1))*((-mu1)*n1m1 + mu2*n1m1 - n1m1*p1 - 2*n2m1*p1 + n1m1*p2 + 
                                           n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))))/
                       sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                       ((-mu1 + mu2 - p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))^2*
                          (mu1^2*mu2*n1m1 - mu1*mu2^2*n1m1 - mu1*mu2*N*p1 + mu2^2*N*p1 + 2*mu1*mu2*n1m1*p1 - 
                             mu2^2*n1m1*p1 + 2*mu1*mu2*n2m1*p1 - mu2*N*p1^2 + mu2*n1m1*p1^2 + 2*mu2*n2m1*p1^2 + 
                             mu1^2*n1m1*p2 - 2*mu1*mu2*n1m1*p2 - 2*mu1*N*p1*p2 + mu2*N*p1*p2 + mu1*n1m1*p1*p2 - 
                             mu2*n1m1*p1*p2 + 2*mu1*n2m1*p1*p2 - mu1*n1m1*p2^2 + mu1*mu2*n1m1*
                             sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) - 
                             mu2*N*p1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                             mu2*n1m1*p1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                             mu1*n1m1*p2*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*tau)/
                       exp((1/2)*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                               4*(mu2*p1 + mu1*(mu2 + p2))))*tau)/(N*p1*(mu2*p1 + mu1*(mu2 + p2))*
                                                                                                     ((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))) - 
                       (exp((1/2)*(-mu1 - mu2 - p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                 4*(mu2*p1 + mu1*(mu2 + p2))))*tau)*(mu1 - mu2 + p1 - p2 + 
                                                                                                       sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))^2*
                          ((-mu1^2)*mu2*n1m1 + mu1*mu2^2*n1m1 + mu1*mu2*N*p1 - mu2^2*N*p1 - 2*mu1*mu2*n1m1*p1 + 
                             mu2^2*n1m1*p1 - 2*mu1*mu2*n2m1*p1 + mu2*N*p1^2 - mu2*n1m1*p1^2 - 2*mu2*n2m1*p1^2 - 
                             mu1^2*n1m1*p2 + 2*mu1*mu2*n1m1*p2 + 2*mu1*N*p1*p2 - mu2*N*p1*p2 - mu1*n1m1*p1*p2 + 
                             mu2*n1m1*p1*p2 - 2*mu1*n2m1*p1*p2 + mu1*n1m1*p2^2 + 
                             mu1*mu2*n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) - 
                             mu2*N*p1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                             mu2*n1m1*p1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                             mu1*n1m1*p2*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*tau)/
                       (N*p1*(mu2*p1 + mu1*(mu2 + p2))*((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))) - 
                       (2*p1*(-1 + (mu1 - mu2 + p1 - p2)/sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                4*(mu2*p1 + mu1*(mu2 + p2))))*
                          (-2 + (mu2*(mu1 + mu2 + p1 + p2 - sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + 
                                                                                                mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                             (1/(N*p1))*(mu1*n1m1 - mu2*n1m1 + n1m1*p1 + 2*n2m1*p1 - n1m1*p2 + 
                                           n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))*tau)/
                       exp((1/2)*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                               4*(mu2*p1 + mu1*(mu2 + p2))))*tau)/sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                                                         4*(mu2*p1 + mu1*(mu2 + p2))) - 
                       (2*exp((1/2)*(-mu1 - mu2 - p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                   4*(mu2*p1 + mu1*(mu2 + p2))))*tau)*p1*(-1 + (-mu1 + mu2 - p1 + p2)/
                                                                                                            sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                          (2 - (mu2*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                  4*(mu2*p1 + mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                             (1/(N*p1))*((-mu1)*n1m1 + mu2*n1m1 - n1m1*p1 - 2*n2m1*p1 + n1m1*p2 + 
                                           n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))*tau)/
                       sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))/
                (2*(4 - (4*mu2*p1)/(mu2*p1 + mu1*(mu2 + p2)) - (4*mu1*p2)/(mu2*p1 + mu1*(mu2 + p2)) - 
                      (2*p1*(-2 + (mu2*(mu1 + mu2 + p1 + p2 - sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                     4*(mu2*p1 + mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                               (1/(N*p1))*(mu1*n1m1 - mu2*n1m1 + n1m1*p1 + 2*n2m1*p1 - n1m1*p2 + 
                                             n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))))/
                      exp((1/2)*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                              4*(mu2*p1 + mu1*(mu2 + p2))))*tau)/sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                                                        4*(mu2*p1 + mu1*(mu2 + p2))) + 
                      ((mu1 - mu2 + p1 - p2 - sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                         (-2 + (mu2*(mu1 + mu2 + p1 + p2 - sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + 
                                                                                               mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                            (1/(N*p1))*(mu1*n1m1 - mu2*n1m1 + n1m1*p1 + 2*n2m1*p1 - n1m1*p2 + 
                                          n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))))/
                      exp((1/2)*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                              4*(mu2*p1 + mu1*(mu2 + p2))))*tau)/sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                                                        4*(mu2*p1 + mu1*(mu2 + p2))) - 
                      (2*exp((1/2)*(-mu1 - mu2 - p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                  4*(mu2*p1 + mu1*(mu2 + p2))))*tau)*p1*
                         (2 - (mu2*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                 4*(mu2*p1 + mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                            (1/(N*p1))*((-mu1)*n1m1 + mu2*n1m1 - n1m1*p1 - 2*n2m1*p1 + n1m1*p2 + 
                                          n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))))/
                      sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))) + 
                      (exp((1/2)*(-mu1 - mu2 - p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                4*(mu2*p1 + mu1*(mu2 + p2))))*tau)*(mu1 - mu2 + p1 - p2 + 
                                                                                                      sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                         (2 - (mu2*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                 4*(mu2*p1 + mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                            (1/(N*p1))*((-mu1)*n1m1 + mu2*n1m1 - n1m1*p1 - 2*n2m1*p1 + n1m1*p2 + 
                                          n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))))/
                      sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))))
  return(gr.mu2f)
}

#######################################################################
#### Graident functions

# Gradient function for Choice Model
gr.choice <- function(par, y){
  grp1 <- gr.p1(par, y)
  grp2 <- gr.p2(par, y)
  grmu1 <- gr.mu1(par, y)
  grmu2 <- gr.mu2(par, y)
  grmat <- matrix(c(grp1, grp2, grmu1, grmu2), ncol = 4)
  return(-colSums(grmat))
}

# Gradient function for p.choice model
gr.p.choice <- function(par, y){
  parf <- c(par[1], par[2], par[3], par[3])
  grp1 <- gr.p1(parf, y)
  grp2 <- gr.p2(parf, y)
  grmu1 <- gr.mu1(parf, y)
  grmu2 <- gr.mu2(parf, y)
  grmat <- matrix(c(grp1, grp2, grmu1, grmu2), ncol = 4)
  grcs <- -colSums(grmat)
  c(grcs[1], grcs[2], sum(grcs[3:4]))
}

# Gradient function for mu.choice model
gr.mu.choice <- function(par, y){
  parf <- c(par[1], par[1], par[2], par[3])
  grp1 <- gr.p1(parf, y)
  grp2 <- gr.p2(parf, y)
  grmu1 <- gr.mu1(parf, y)
  grmu2 <- gr.mu2(parf, y)
  grmat <- matrix(c(grp1, grp2, grmu1, grmu2), ncol = 4)
  grcs <- -colSums(grmat)
  c(sum(grcs[1:2]), grcs[3], grcs[4])
}

# Gradient function for Fixed model
gr.fixed <- function(par, y){
  parf <- c(par[1], par[1], par[2], par[2])
  grp1 <- gr.p1(parf, y)
  grp2 <- gr.p2(parf, y)
  grmu1 <- gr.mu1(parf, y)
  grmu2 <- gr.mu2(parf, y)
  grmat <- matrix(c(grp1, grp2, grmu1, grmu2), ncol = 4)
  grcs <- -colSums(grmat)
  c(sum(grcs[1:2]), sum(grcs[3:4]))
}

#### Verify gradient function with numerical approximation
# Should be commented out to "source" the file
# Need to define pars, a vector of true parameter values,
# yt, a dataset with t, n1, n2, and n3 columns
# and NLL.choice(), the likelihood function

# # # Data set
# yt <- structure(list(t = 1:20, n1 = c(12L, 15L, 17L, 15L, 16L, 17L,
#                                       19L, 19L, 16L, 18L, 20L, 19L, 19L, 19L, 18L, 19L, 19L, 18L, 20L,
#                                       19L), n2 = c(2L, 1L, 3L, 2L, 2L, 1L, 0L, 0L, 0L, 0L, 0L, 0L,
#                                                    0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L), n3 = c(2L, 0L, 0L, 3L, 1L, 1L,
#                                                                                            1L, 1L, 3L, 2L, 0L, 1L, 0L, 1L, 1L, 1L, 1L, 2L, 0L, 1L), `NA` = c(4L,
#                                                                                                                                                              4L, 0L, 0L, 1L, 1L, 0L, 0L, 1L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L,
#                                                                                                                                                              0L, 0L, 0L), N = c(20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20,
#                                                                                                                                                                                 20, 20, 20, 20, 20, 20, 20, 20, 20)), .Names = c("t", "n1", "n2",
#                                                                                                                                                                                                                                  "n3", NA, "N"), row.names = c(NA, -20L), class = "data.frame")
# # Paramter value vector
# p1t = 0.368; p2t = p1t; mu1t = 0.167; mu2t = 0.0001 # Old param estimates for Brown-TBW 2010 data
# pars <- c(p1t, p2t, mu1t, mu2t)
# 
# pars <- rep(0.01, 4)
# yt <- dat[,c("t", "n1", "n2", "n3")]
# # # Compare gradient equation to numerical approximation
# library(numDeriv)
# exact.choice <- gr.choice(pars, yt) # Gradient equation
# numd.choice <- grad(function(u) NLL.choice(u, yt), pars) # Numerical approximated gradient
# rbind(exact.choice, numd.choice)
