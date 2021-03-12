### Neutrophils self-limit swarming to contain bacterial growth in vivo

library(ggplot2)

### fig. S1G
names(data)  
# output: "Position.X", "Position.Y","Position.Z","Unit","Category","Collection","Time","Parent" (= ID of track),"ID" (= Id of spots),"X" (=WT or KO)
# -Y gives direction to attractor well

getDerivatives <- function(data){
  dataGetDer <- NULL
  for(track in levels(data$Parent)){
    objectsInTrack <- data[data$Parent==track,]
    objectsInTrack <- objectsInTrack[order(objectsInTrack[,"Time"]),]
    
    dt <- diff(objectsInTrack[,"Time"])
    dx <- diff(objectsInTrack[,"Position.X"])
    dy <- diff(objectsInTrack[,"Position.Y"])
    dxy <- sqrt(dx^2+dy^2)
    cum_dxy <- cumsum(dxy)
    vx <- dx/dt
    vy <- dy/dt
    v <- sqrt(vx^2+vy^2)
    
    objectsInTrack <- objectsInTrack[-1,]    # drops first point (incomplete derivative information)
    chemIdx <- -(dy / dxy) #"velocity angle" = acos(chemIdx)*180/pi
    chemIdx[is.nan(chemIdx)] <- 0
    objectsInTrack <- cbind(objectsInTrack,dx,dy,dxy,cum_dxy,vx,vy,v,chemIdx)
    dataGetDer <- rbind(dataGetDer,objectsInTrack)
  }
  return(dataGetDer)
}

# for: WT is desensitized after 90min
data_91 <- data[data$X=="WT",]                                  
data_startingpoints_91 <- data_91[data_91$Time==91,c(1,2,8)]       
colnames(data_startingpoints_91) <- c("StartPos.X", "StartPos.Y","Parent")    
data_91 <- merge(data_startingpoints_91, data_91, by="Parent")
data_91$DistanceXStartpoint <- (data_91$Position.X - data_91$StartPos.X)        
data_91$DistanceYStartpoint <- (data_91$Position.Y - data_91$StartPos.Y)        
data_91$DistanceXY <- sqrt((data_91$DistanceYStartpoint)^2 + (data_91$DistanceXStartpoint)^2)    
data_91 <- data_91[data_91$Time >90,]
rm(data_startingpoints_91)
dataGetDer_91 <- getDerivatives(data=data_91)

theme_set(theme_bw())
ggplot(data=dataGetDer_91,aes(x=Time,y=DistanceYStartpoint)) + scale_y_reverse() +  
  geom_line(aes(colour=chemIdx,group=Parent),size=0.5 ) +
  scale_colour_gradient2(name="Chemotactic\nIndex",low="purple4", mid="magenta", high="orange",  space="rgb",limits=c(-1,1))+
  facet_wrap(~X)



### calculations in Fig. 1C, Fig. 1D, and fig. S1G
# "start" and "end" define timeframe of interest
speed <- (dataGetDer[dataGetDer$Time == end,"cum_dxy"] - dataGetDer[dataGetDer$Time == start,"cum_dxy"] ) / (end-start)
y-straightness <- -1* (dataGetDer[dataGetDer$Time == end,"Position.Y"] - dataGetDer[dataGetDer$Time == start,"Position.Y"]) / (dataGetDer[dataGetDer$Time == end,"cum_dxy"] - dataGetDer[dataGetDer$Time == start,"cum_dxy"] )



### fig. S3D-F
names(data)  
# output: "Position.X", "Position.Y","Position.Z","Unit","Category","Collection","Time","Parent" (= ID of track),"ID" (= Id of spots),"X" (=WT or KO)data <- data[data$X=="WT",]

data_startingpoints <- data[data$Time==1,c(1,2,8)]       
colnames(data_startingpoints) <- c("StartPos.X", "StartPos.Y","TrackID")    
data <- merge(data_startingpoints, data, by="TrackID")
data$DistanceXStartpoint <- (data$Position.X - data$StartPos.X)       
data$DistanceYStartpoint <- (data$Position.Y - data$StartPos.Y)        
data <- data [order(data[,"Time"]),]                                   
rm(data_startingpoints)

theme_set(theme_bw())
ggplot(data,aes(x=DistanceXStartpoint,y=DistanceYStartpoint)) +
  geom_path(aes(color=TrackID, group=TrackID),size=1 ) +
  xlim(-200,200) + ylim (-200,200) + 
  theme(legend.position = "none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))



