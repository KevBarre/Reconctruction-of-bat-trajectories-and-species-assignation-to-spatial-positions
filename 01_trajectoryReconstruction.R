# This script allows to gather spatial bat positions belonging to a same trajectory within an 
# identifier unique to this trajectory. Several metrics such as speed or height of flight,
# and the angle to an object (a light in this case) are then computed.
#############################################################################################
# Written by Yves Bas, Julie Pauwels, Charlotte Roemer and Kévin Barré
#############################################################################################

rm(list=ls())

library(MASS) # V7.3-50
library(data.table) # V1.11.4
library(gdata) # V2.18.0
library(beepr) # V1.3
library(stringr) # V1.3.1

# working directory
setwd("./")

# if all your trajectory files are already concatenated please disable L14-31 and activate L32
# folder containing the .csv of the positions recorded
# dir.pos = "C:/Users/barre/Documents/Postdoc_Chirolum/Papiers/PaysBas_trajecto/DATA/FILES_XLS"
# 
# # list  of .csv files names
# filenames.pos = list.files(dir.pos, pattern=".xls", full.names=T, recursive=F)
# # lapply(filenames.pos, function(f) {
# #   df = read.xlsx(f, sheet=1)
# #   write.csv(df, gsub("xls", "csv", f), row.names=FALSE)
# # })
# # creates a list containing all the tables of positions
# tablelist.pos = list()
# for (i in 1:length(filenames.pos))
# {
#   tablelist.pos[[i]] = read.xls(filenames.pos[[i]], skip=3, perl = "C:/Perl64/bin/perl.exe")
#   tablelist.pos[[i]] = cbind(fichier = filenames.pos[i], tablelist.pos[[i]])
# }

# concatenates all the tables of the list to a unique table
# table.pos = as.data.frame(rbindlist(tablelist.pos))

table.pos = fread("./01_dataSample.csv")
# This file must have following columns at least:
#   - "fichier" = trajectory file positions come from (i.e. if several file were recorded a same night and/are if several night were recorded)
#   - "temps..s." = time of the bat call since the recording start
#   - "x..m." = x position in meters
#   - "X.dx..m..1" = imprecision of the x position in meters
#   - "y..m." = y position in meters
#   - "X.dy..m..1" = imprecision of the y position in meters
#   - "y..m." = y position in meters
#   - "X.dy..m..1" = imprecision of the y position in meters
#   - "z..m." = z position in meters
#   - "X.dz..m..1" = imprecision of the z position in meters
#   - "fsigmax..Hz." = energy peak of calls
#   - "ID" = identifier of the sampling site in a given night

# change numbers within the table as numeric
for (i in 2:18)
{
  table.pos[,i] = as.numeric(as.character(table.pos[, i]))
}

### FUNCTIONS ##########################################################################################################

# SeqTraj : separates in the table of positions in blocks. A new block is defined when two position are more 
# than 2 seconds appart.
# x = vector of the time of recording of each position
# y = vector of the file name for each position
SeqTraj = function(x, y)
{
  IdTraj = rep(1, length(x)) # vector that will contain blocks identifiers
  for (i in 2:length(x))
  {
    if (y[i] == y[i-1]) # if the position is in the same file as the previous (i.e. same night of recording)
    {
      if (x[i] - x[i-1] > 2) # if the time difference is above 2 seconds
      {
        IdTraj[i] = IdTraj[i-1] + 1 # this position's identifier is the previous position's identifiers + 1
      }
      else # if the time difference is under 2 seconds
      {
        IdTraj[i] = IdTraj[i-1] # this position's identifier is the same as the previous position's identifiers
      }
    }
    else # if the position is in a different file than the previous (i.e. different night of recording)
    {
      IdTraj[i] = IdTraj[i-1] + 1# this position's identifier is the previous position's identifiers + 1
    }
  }
  IdTraj
}

# Dist : compute the tridimentional distance between 2 positions
Dist = function(p1, p2)
{
  ((p2$x..m.-p1$x..m.)^2 + (p2$y..m.-p1$y..m.)^2 + (p2$z..m.-p1$z..m.)^2)^0.5
}

# DAC = distance to the "center" of the recording antenna (0,0)
DAC = function(p)
{
  ((0 - p$x..m.)^2 + (0 - p$y..m.)^2 + (0 - p$z..m.)^2)^0.5
}

# Vitesse : speed between two positions
Vitesse = function(p1,p2)
{
  if (p2$temps..s. == p1$temps..s.) # if time of recording is identical
  {
    999
  }
  else
  { # distance between the points / time difference of recording
    ((p2$x..m. - p1$x..m.)^2 + (p2$y..m. - p1$y..m.)^2 + (p2$z..m. - p1$z..m.)^2)^0.5 / (p2$temps..s. - p1$temps..s.)
  }
}

# Vhor : horizontal speed, i.e. speed in the (x,y) plan
Vhor=function(p1, p2)
{ # distance between the points considering only x and y / time difference of recording
  ((p2$x..m. - p1$x..m.)^2 + (p2$y..m. - p1$y..m.)^2)^0.5 / (p2$temps..s. - p1$temps..s.)
}

# CoefVar : coefficient of variation
CoefVar=function(x)
{
  sd(x)/mean(x) 
}

### TABLE PREPARATION TO TRAJECTORIES RECONSTRUCTION ##################################################################

# Compute the block identifier of each position
IdTraj = SeqTraj(table.pos$temps..s., table.pos$fichier)

# Measure of positions imprecision
# Squared root of the sum of squared imprecisions for the 3 dimensions
Imp = ((table.pos$X.dx..m..1)^2 + (table.pos$X.dy..m..1)^2 + (table.pos$X.dz..m..1)^2)^0.5 
table.pos = cbind(table.pos, IdTraj, Imp)
# Inter-pulse duration of consecutive positions
D = rep(999, nrow(table.pos))
for (i in 2:nrow(table.pos))
{
  if (IdTraj[i] == IdTraj[i-1]) # if the position is in the same block than the previous one
  {
    D[i] = table.pos$temps..s.[i]-table.pos$temps..s.[i-1] # Compute the time between positions
  }
}
table.pos = cbind(table.pos, D)
table.pos = subset(table.pos,Imp < 1) # Remove positions for which imprecision is > 1 meter

# # Initialize a speed vector
V = rep(999, nrow(table.pos))
# Compute the speed between adjacent positions within a block
for (i in 2:nrow(table.pos))
{
  if (IdTraj[i] == IdTraj[i-1]) # if the position is in the same block as the previous one
  {
    V[i] = Vitesse(table.pos[i-1, ], table.pos[i, ]) # Compute the speed
  }
}
# Plot the speed versus block identifiers
plot(V, table.pos$IdTraj, xlim=c(0, 150))

# Initialize a speed coefficient of variation vector
CvV = rep(999, nrow(table.pos))
# Coefficient of variation for 3 consecutive speeds (4 points). Stored at the third position
for (i in 3:(nrow(table.pos)-1))
{
  if (IdTraj[i+1]==IdTraj[i-2]) # if the fourth position is in the same block as the first one
  {
    CvV[i] = CoefVar(V[(i-1):(i+1)])  # Compute de speed coefficient of variation
  }
}

# Plot the coefficient of variation versus the block identifiers
plot(CvV, table.pos$IdTraj, xlim=c(0, 4))

# 
# # Idem for horizontal speed
# VH = rep(999,nrow(table.pos))
# for (i in 2:nrow(table.pos))
# {
#   if (IdTraj[i] == IdTraj[i-1]) 
#   {
#   VH[i] = Vhor(table.pos[i-1,],table.pos[i,])
#   }
# }
# 
# CvVh = rep(999,nrow(table.pos))
# for (i in 3:(nrow(table.pos)-1))
# {
#   if (IdTraj[i+1] == IdTraj[i-2]) 
#   {
#   CvVh[i] = CoefVar(VH[(i-1):(i+1)])
#   }
# }

# Compute flight height, distance to the light and angle accounting for antena height
# A negative angle means the position is under the light
date = c("20180710", "20180711", "20180712", "20180713", "20180714", "20180716", "20180717", "20180718", "20180719",
         "20180720", "20180721", "20180722")
antena = c(0.95, 0.90, 0.98, 0.99, 0.82, 1, 0.95, 0.90, 0.94, 0.91, 0.95, 0.98) # antenna height in cm
height = rep(999, nrow(table.pos))
distLight = rep(999, nrow(table.pos))
angle = rep(999, nrow(table.pos))
pb <- txtProgressBar(min = 0, max = nrow(table.pos), style = 3)
for (i in date) 
  {
  print(i)
  for (j in 2:nrow(table.pos)) 
    {
    Sys.sleep(0.1)
    setTxtProgressBar(pb, j)
    if (str_detect(table.pos$ID, i)) 
      {
      height[j] = table.pos$z..m.[j] + antena[which(date==i)]
      distLight[j] = ((table.pos$x..m.[j]+4)^2 + (table.pos$y..m.[j]+4)^2 + (height[j]-5)^2)^0.5
      angle[j] = (asin((height[j]-5)/distLight[j])*180)/3.141592654
    }
  }
}

# Add the block identifier to the table of positions (V & CVv are not added because they are calculated within the TriTraj function)
table.pos2=cbind(table.pos, V, CvV, height, distLight, angle)
beep(8)
write.csv(table.pos2, "./TriTrajBrut.csv")

### TRAJECTORIES RECONSTRUCTION ########################################################################################

# TriTraj : reconstruct one trajectory per block (if possible)
# tab = table of positions
# r = round 
TriTraj = function(tab, r)
{
  TFTot = tab[0,] # Trajectoires Filtrées Total : table containing blocks of positions constituting a trajectory (blocks
                #identified with IdTraj + round)
  PositOrph = tab[0,] # Positions Orphelines : table containing all the positions not used in any trajectory of TFTot
                    # (trajectory not reconstituted yet or useless position)
  
  for (j in 1:max(tab$IdTraj)) # j : block identifier
  {
    print(paste0("round ", r, " block ", j, "/4273"))
    subtest = subset(tab,tab$IdTraj == j) # subset of position table for identifier j
    Fmed_temp = quantile(subtest$fsigmax..Hz., 0.5)
    # Positions for which the energy peak = median of the block's energy peaks +/- 5 kH
    TrajFiltr0 = subset(subtest, abs(subtest$fsigmax..Hz. - Fmed_temp) < 5000) 
    
    # Positions for which the energy peak is outside this range
    HorsFiltr = subset(subtest,abs(subtest$fsigmax..Hz. - Fmed_temp) >= 5000) 
    
    if (nrow(TrajFiltr0) > 4) # if there are more than 4 positions (else useless)
    {
      # Compute the speed
      # V = rep(999, nrow(TrajFiltr0)) 
      # for (k in 2:nrow(TrajFiltr0)) {V[k] = Vitesse(TrajFiltr0[k-1, ], TrajFiltr0[k, ])}
      # # Compute the coefficient of variation of the speed
      # CvV = rep(999, nrow(TrajFiltr0))
      # for (k in 3:(nrow(TrajFiltr0)-1)){ CvV[k] = CoefVar(V[(k-1):(k+1)])}
      
      # table containing the 3 consecutive lines that generate the smallest coefficient of variation
      TrajFiltr = TrajFiltr0[max(1, (which.min(TrajFiltr0$CvV)-2)):(which.min(TrajFiltr0$CvV)), ]
      # TrajFiltr = TrajFiltr0[(which.min(CvV)-1):(which.min(CvV)+1), ]
      
      if (which.min(TrajFiltr0$CvV)>3) # if the 3 lines selected are not the first 3 lines of TrajFiltr0...
      # if (which.min(CvV)>2) 
          
      {
        for (i in (which.min(TrajFiltr0$CvV)-3):1) # ...go back up the block one line at the time
        # for (i in (which.min(CvV)-2):1)
        {
          d = Dist(TrajFiltr0[i,],TrajFiltr[1,]) # compute the distance between the first line of TrajFiltr0 and the line above
          dt = (TrajFiltr$temps..s.[1]-TrajFiltr0$temps..s.[i]) # compute the time difference between the first line of TrajFiltr0 and the line above
          
          if (( # check several conditions before adding the position to TrajFiltr:
            Vitesse(TrajFiltr0[i,],TrajFiltr[1,]) < Vitesse(TrajFiltr[1,],TrajFiltr[2,]) + Vitesse(TrajFiltr[2,],TrajFiltr[3,]) -dt-d) 
            # 1) the speed of the line to add is smaller than the sum of the sum of the speeds of the two following lines
            # minus penalties if the distance or time difference is important
            & (dt < 2 - (quantile(subtest$fsigmax..Hz., 0.5)>30)) # 2) there are no more than 2 seconds between the position to add and the following one
            & (d/dt < 30) ) # 3) the speed is below 30 m/s
          {
            TrajFiltr = rbind(TrajFiltr0[i, ], TrajFiltr) # add the position to TrajFiltr
          }else{
            if (d > 3) {HorsFiltr = rbind(TrajFiltr0[i, ], HorsFiltr) } 
            # For positions that are do not check the conditions to be part of the trajectory (TrajFiltr), keep those that are at least 3 meters away from the first of TrajFiltr0
            # Hypothesis : the positions closer than 3 meters have great chances to be echoes of the trajectory in Trajfiltr. Removing them avoid reconstructing "ghost" trajectories.
          }
        }
      }
      
      if (which.min(TrajFiltr0$CvV)!=nrow(TrajFiltr0)) # if the 3 lines selected are not the last 3 lines of TrajFiltr0...
      # if ((which.min(CvV)+1)!=nrow(TrajFiltr0))
      {
        
        for (h in (which.min(TrajFiltr0$CvV)+1):nrow(TrajFiltr0)) # ...go down the block one line at the time
        # for (h in (which.min(CvV)+2):nrow(TrajFiltr0))
        {
          d = Dist(TrajFiltr[nrow(TrajFiltr),],TrajFiltr0[h,])
          dt = (TrajFiltr0$temps..s.[h]-TrajFiltr$temps..s.[nrow(TrajFiltr)])
          
          if (
            ((Vitesse(TrajFiltr[nrow(TrajFiltr),],TrajFiltr0[h,]))
             <
             (if (nrow(TrajFiltr)>2) {Vitesse(TrajFiltr[nrow(TrajFiltr)-1,],TrajFiltr[nrow(TrajFiltr),])}
              else{30})
             + (if (nrow(TrajFiltr)>2) {Vitesse(TrajFiltr[nrow(TrajFiltr)-2,],TrajFiltr[nrow(TrajFiltr)-1,])}
                else{30})
            ) & (dt < 2- (quantile(subtest$fsigmax..Hz., 0.5)>30))
            & (d/dt<30)
          )
          {
            TrajFiltr=rbind(TrajFiltr,TrajFiltr0[h,])
          }else{
            if (d > 3) {HorsFiltr = rbind(HorsFiltr, TrajFiltr0[h, ])}
          }
        }
      }
      
      plot(TrajFiltr$x..m.,TrajFiltr$y..m.,xlim=c(-25,25),ylim=c(-15,15),type="l") # 2D plot of the trajectory reconstructed
      
      TFTot=rbind(TFTot,TrajFiltr) # store the trajectory
      
      PositOrph=rbind(PositOrph,HorsFiltr) # store the unused positions for following rounds
      
    }
  }
  
  if (nrow(TFTot)>0)
  {
    Pool=cbind(TFTot,round=r) # All trajectories of TFTot are assigned a round number
    
    Pool2=cbind(PositOrph,round=r+1) # All unused positions are assigned the round number +1 (to be analysed in the next round)
    
    rbind(Pool,Pool2) # Bind trajectories already reconstructed and positions that could participate in one
                      # solitary positions (TrajFiltr0<3) and echo positions are eliminated
  } else
    {
    tabtemp=cbind(tab,round=r)
    tabtemp[0,]
  }
}

TriTrajT=TriTraj(table.pos2,1)
TriTrajT=TriTrajT[order(TriTrajT$temps..s.),]

ATrier=subset(TriTrajT,TriTrajT$round==max(TriTrajT$round))
AGarder=subset(TriTrajT,TriTrajT$round<max(TriTrajT$round))

# All unused positions are re-launched, then the next round re-launches those not used by the previous round, etc...
for (j in 3:10){
while (nrow(ATrier)>5) {
  j=j+1
  print(j)
  
  Vt=rep(999,nrow(ATrier))
  
  for (i in 2:nrow(ATrier))
  {
    
    if (ATrier$IdTraj[i]==ATrier$IdTraj[i-1]) 
    {
      Vt[i]=Vitesse(ATrier[i-1,],ATrier[i,])
    }
  }
  
  CvVt=rep(999,nrow(ATrier))
  
  for (i in 3:(nrow(ATrier)-1))
  {
    if (ATrier$IdTraj[i+1]==ATrier$IdTraj[i-2]) 
    {
      
      CvVt[i]=CoefVar(Vt[(i-1):(i+1)])
    }
  }
  
  ATrier$CvV=CvVt
  TriTrajSub=TriTraj(ATrier[,1:(ncol(ATrier)-1)],j)
  TriTrajSub=TriTrajSub[order(TriTrajSub$temps..s.),]
  AGarder=rbind(AGarder,subset(TriTrajSub,TriTrajSub$round==j))
  ATrier=subset(TriTrajSub,TriTrajSub$round==j+1)
  
}

AGarder=AGarder[order(AGarder$IdTraj,AGarder$round,AGarder$temps..s.),]

# Updating of speed flight
Vfinal=rep(999,nrow(AGarder))

for (i in 2:nrow(AGarder))
{
   if((AGarder$IdTraj[(i-1)]==AGarder$IdTraj[i])&(AGarder$round[(i-1)]==AGarder$round[i]))
   {
    Vfinal[i]=Vitesse(AGarder[(i-1),],AGarder[i,])
    
  }
}

Vmed=aggregate(Vfinal,by=list(AGarder$IdTraj,AGarder$round),FUN=function(x) quantile(x,0.5))
head(AGarder)

AGarderFinal<-cbind(AGarder,Vfinal)

# Saving
write.table(AGarderFinal,"./TriTrajAGarder.csv",sep=";",row.names=F)

}

# Adding all individual positions
table.pos2 <- fread("./TriTrajBrut.csv")
AGarderFinal <- fread("./TriTrajAGarder.csv")

Agarder_new = AGarderFinal
table.pos3 = table.pos2
round =rep(NA,nrow(table.pos3))
table.pos3 = cbind(table.pos3, round)
table.pos3$Vfinal = table.pos3$V

table.pos3$IdTraj = paste("i", table.pos3$IdTraj, "_",table.pos3$V1, sep = "")

pb <- txtProgressBar(min = 0, max = nrow(table.pos2), style = 3)

for (i in 1:nrow(table.pos3)) {
  Sys.sleep(0.1)
  setTxtProgressBar(pb, i)
  if (sum(which(AGarderFinal$V1==table.pos3$V1[i])) == 0) {
    Agarder_new = rbind(Agarder_new, table.pos3[i])
  }
}

write.table(Agarder_new,"C:/Users/barre/Documents/Postdoc_Chirolum/Papiers/PaysBas_trajecto/DATA/TriTrajAGarderFull.csv",sep=";",row.names=F)

