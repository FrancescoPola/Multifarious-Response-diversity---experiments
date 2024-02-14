##
# ad hoc function for the exploration analysis assessing the effects of diversity in z

plot_RD <- function(dir_1, dir_2, dir_3, rdiv_1, rdiv_2, rdiv_3) {

  rmax<-max(rdiv_1$rdiv,rdiv_2$rdiv,rdiv_3$rdiv); rmin<-1
  dmax<-max(rdiv_1$sign,rdiv_2$sign,rdiv_3$sign); dmin<-0
  
  Fig_1_comm1 <-dir_1 %>%
    ggplot(aes(x=time, y=dir_deriv, col = sp)) +
    #theme_classic(base_size = 14) + 
    labs(x = "time",y = "Directional derivative",tag = "a)") + 
    #geom_point(size=0.5) +
    geom_line() + theme_bw() + geom_hline(yintercept = 0, linetype= "dashed")
  
  Fig_1_comm2 <-dir_2 %>%
    ggplot(aes(x=time, y=dir_deriv, col = sp)) +
    #theme_classic(base_size = 14) + 
    labs(x = "time",y = "Directional derivative",tag = "b)") + 
    #geom_point(size=0.5) +
    geom_line() + theme_bw() + geom_hline(yintercept = 0, linetype= "dashed")
  
  Fig_1_comm3 <-dir_3 %>%
    ggplot(aes(x=time, y=dir_deriv, col = sp)) +
    #theme_classic(base_size = 14) + 
    labs(x = "time",y = "Directional derivative",tag = "c)") + 
    #geom_point(size=0.5) +
    geom_line() + theme_bw() + geom_hline(yintercept = 0, linetype= "dashed") 
  
  
  Fig2_comm1 <- rdiv_1 %>% 
    ggplot(mapping = aes(x = time,y = rdiv)) +
    theme_bw(base_size = 12) + 
    labs(x = NULL,y = "dissimilarity (derivatives)",tag = "d)") + 
    geom_hline(yintercept = rdiv_1$Med,lty=2) +
    geom_line() + 
    geom_richtext(x = 7,
                  mapping = aes(y = rmax),
                  size=4.5,
                  label.color = NA,
                  label = paste("RDiv^Median =",paste0(round(rdiv_1$Med,digits = 3))))+
    lims(y = c(rmin,rmax))
  
  Fig2_comm2 <- rdiv_2 %>% 
    ggplot(mapping = aes(x = time,y = rdiv)) +
    theme_bw(base_size = 12) + 
    labs(x = NULL,y = "dissimilarity (derivatives)",tag = "e)") + 
    geom_hline(yintercept = rdiv_2$Med,lty=2) +
    geom_line() + 
    geom_richtext(x = 7,
                  mapping = aes(y = rmax),
                  size=4.5,
                  label.color = NA,
                  label = paste("RDiv^Median =",paste0(round(rdiv_2$Med,digits = 3))))+
    lims(y = c(rmin,rmax))
  
  Fig2_comm3 <- rdiv_3 %>% 
    ggplot(mapping = aes(x = time,y = rdiv)) +
    theme_bw(base_size = 12) + 
    labs(x = NULL,y = "dissimilarity (derivatives)",tag = "f)") + 
    geom_hline(yintercept = rdiv_3$Med,lty=2) +
    geom_line() + 
    geom_richtext(x = 7,
                  mapping = aes(y = rmax),
                  size=4.5,
                  label.color = NA,
                  label = paste("RDiv^Median =",paste0(round(rdiv_3$Med,digits = 3))))+
    lims(y = c(rmin,rmax))
  
  
  Fig3_comm1 <- rdiv_1 %>% 
    ggplot(mapping = aes(x = time,y = sign)) +
    theme_bw(base_size = 12)+
    labs(x = "time",y = "Divergence", tag = "g)") + 
    geom_line() +
    geom_hline(yintercept = rdiv_1$Med_sign,lty=2) +
    geom_line() + 
    geom_richtext(x = 7,
                  mapping = aes(y = dmax),
                  size=4.5,
                  label.color = NA,
                  label = paste("RDiv^Median =",paste0(round(rdiv_1$Med_sign,digits = 2))))+
    lims(y = c(dmin,dmax))
  
  Fig3_comm2 <- rdiv_2 %>% 
    ggplot(mapping = aes(x = time,y = sign)) +
    theme_bw(base_size = 12)+
    labs(x = "time",y = "Divergence", tag = "h)") + 
    geom_line() +
    geom_hline(yintercept = rdiv_2$Med_sign,lty=2) +
    geom_line() + 
    geom_richtext(x = 7,
                  mapping = aes(y = dmax),
                  size=4.5,
                  label.color = NA,
                  label = paste("RDiv^Median =",paste0(round(rdiv_2$Med_sign,digits = 2))))+
    lims(y = c(dmin,dmax))
  
  Fig3_comm3 <- rdiv_3 %>% 
    ggplot(mapping = aes(x = time,y = sign)) +
    theme_bw(base_size = 12)+
    labs(x = "time",y = "Divergence", tag = "i)") + 
    geom_line()+
    geom_hline(yintercept = rdiv_3$Med_sign,lty=2) +
    geom_line() + 
    geom_richtext(x = 7,
                  mapping = aes(y = dmax),
                  size=4.5,
                  label.color = NA,
                  label = paste("RDiv^Median =",paste0(round(rdiv_3$Med_sign,digits = 2))))+
    lims(y = c(dmin,dmax))
  
  
  
  (Fig_1_comm1 + Fig_1_comm2 + Fig_1_comm3 )/ (Fig2_comm1 + Fig2_comm2 + Fig2_comm3 )/ (Fig3_comm1 + Fig3_comm2 + Fig3_comm3 )
  
}

