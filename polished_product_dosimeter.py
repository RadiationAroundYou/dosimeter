#Copyright Steve Lloyd, Tom Stevenson and the Langton Star Centre




import glob
import sys
import re
import os
import struct
import string
import math
import cgi
import cgitb

exposure_list = []
metadata_list = []
dose_list = []

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

#Chip Parameters:

chip_volume = (0.0002) * (0.0002) #Volume of the chip (in m^3), area of 2cm^2 by 200 micrometers deep
chip_density = 2330                    #Density of the chip (in kg/(m^3)) (Taking that of silicon for now)
chip_mass = chip_volume * chip_density

metadata_position = 153                  # This is the number of bytes into the metadata file the acq time line is on
start_time_position = 1380                  #Ditto for the start time (UNIX time)
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

#4

def read_in_dir(path):

    
    listing = os.listdir(path)

    
    for f in listing:
        file_object = open(path + '/' + f)
        first_chars = file_object.read(2)
        file_object.close()
        if first_chars == 'A0' :        #Finds if it's metadata
            process_metadata(path + '/' + f)

            
    for f in listing:
        file_object = open(path + '/' + f)
        first_chars = file_object.read(1)
        file_object.close()
        if first_chars == '0' :    #Finds if it's a pixelman output data file
            for item in metadata_list:
                 if item[0][:-4] == (f): #Matches the exposures to the data files
                    exposure_time = item[1] #Takes the exposure from the list
                    start_time = item[2]    #Takes the  start time from the list
                    read_pix_clusters_dose(path + '/' + f, exposure_time)
                    dose_list[-1].append(start_time)
   


    total_dose_rate = 0

    for element in dose_list:
        print '=============================='
        print element[2]
        print 'Filename:', element[0]
        print 'Aquisition time:', element[1], 's'
        print 'Start time (UNIX):', element[6]
        print element[2] * (10 ** 12 ), 'pSv per second'
        print ( element[2] * 3600 ) * (10 ** 9), 'nSv per hour'
        print (element[2] * 3600 * 24 * 365 ) * (10 ** 3), 'mSv per year'
        print element[3] * (10 ** 12), 'pSv per second (Alpha)'
        print element[4] * (10 ** 12), 'pSv per second (Beta)'
        print element[5] * (10 ** 12), 'pSv per second (Gamma)'
        total_dose_rate += (element[2])

    average_dose_rate = total_dose_rate / len(dose_list) 
    print '================================='
    print average_dose_rate * 3600 * (10 ** 9), 'nSv per hour average'   
        
    write_output(dose_list)    #Prints output in CSV format: quicker than printing to standard output
    
                     
    
def process_metadata(pathname):

        filename = pathname.split('/')[-1]


        acq_time = get_acq_time(pathname)
        start_time = get_start_time(pathname)
        metadata_list.append([filename, acq_time, start_time]) 


def get_acq_time(pathname):
    f = open(pathname, 'r')
    acq_flag = 3
    for lines in f:
        acq_flag += 1
        if acq_flag == 2:
            acq_time = float(lines.strip()) #Removes newline char
            return acq_time
        if (string.find(lines, 'Acq time') != -1):
             acq_flag = 0     #line after next contains the value of the acq time
    f.close()    

def get_start_time(pathname):
    f = open(pathname, 'r')
    time_flag = 3         
    for lines in f:
        time_flag += 1
        if time_flag == 2:
            start_time = float(lines.strip())
            return start_time
        if string.find(lines, 'Start time') != -1:
             time_flag = 0
    f.close()         
             
            
            



def read_pix_clusters_dose(pathname, exposure_time):
    
    filename = pathname.split('/')[-1]

    hits = get_hits(pathname)
    global isBackground
    if (total_sum(hits) / exposure_time) > 100: 
        isBackground = False    #I think this works....
    else:
        isBackground = True


    clusters = get_clusters(hits)

    (alpha_doserate, beta_doserate, gamma_doserate, total_doserate) = get_doserates(clusters, exposure_time)


    dose_list.append([filename, exposure_time, total_doserate, alpha_doserate, beta_doserate, gamma_doserate])

def get_doserates(clusters, exposure_time):
    
    total_dose_rate = 0
    alpha_dose_rate = 0
    beta_dose_rate = 0
    gamma_dose_rate = 0
    
    for cluster in clusters:
        if cluster.type == 'Alpha':
            alpha_dose_rate += (cluster.equivalentenergy / chip_mass) / exposure_time
        if cluster.type == ('Beta' or 'Maaybe Alpha'):
            beta_dose_rate += (cluster.equivalentenergy / chip_mass) / exposure_time
        if cluster.type == 'Gamma':
            gamma_dose_rate += (cluster.equivalentenergy / chip_mass) / exposure_time

        total_dose_rate += (cluster.equivalentenergy / chip_mass) / exposure_time    
    
    return (alpha_dose_rate, beta_dose_rate, gamma_dose_rate, total_dose_rate)


def write_output(dose_list):
    output = open('csv_output.csv', 'w')
    output.write('Filename,Acquisition Time,Total Dose Rate,Alpha Dose Rate,Gamma Dose Rate,Beta Dose Rate,Gamma Dose Rate,Start Time (UNIX Time)\n')

    output_list = dose_list[:]
    for element in output_list:
        for item in element:
            item = str(item)
            print item    #Changes all into strings so the .join method can be used

    output_list = []
    for element in dose_list:
        output_sub_list = []
        for item in element:
            output_sub_list.append(str(item))
        output_list.append(output_sub_list)
        

    for element in output_list:
        print element
        output.write(",".join(element))
        output.write('\n')
    output.close()    
            
    
        
        
                            
        
#-----------------------------------------------------------------------------
def get_hits(pathname):

    # Reads the input file and returns a list of all hits.
    
    f = open(pathname, 'r')
    text = f.read()
    f.close()
    
    hits = []

                

    rows = text.splitlines()
    i = 0
    for row in rows:
#       print i, row
        vals = row.split(' ')
        j = 0
        adjacent_neighbours = 0
        count = hits.count
        for val in vals:
            if val != '0':
                hit = Hit(i, j, val, count, adjacent_neighbours )
                hits.append(hit)
#               print i,j, val
            j += 1
        i += 1
    
    return hits

#-----------------------------------------------------------------------------

def get_neighbours(hits):
    neighbours_list = []
    for hit in hits:
        neighbours = 0  #Iterates through hits, finds number of neighbours and returns a list of numbers of neighbours
        for nhit in hits:
            if abs(hit.i-nhit.i) <= 1 and abs(hit.j-nhit.j) <= 1: 
                neighbours += 1
        neighbours -= 1
        neighbours_list += str(neighbours)
    return neighbours_list


#-----------------------------------------------------------------------------

def get_adjacent_neighbours(hits):
    neighbours_list = []
    for hit in hits:
        neighbours = 0
        for nhit in hits:
            if abs(hit.i-nhit.i) <= 1 and abs(hit.j-nhit.j) <= 1 and (abs(hit.i-nhit.i) + abs(hit.j-nhit.j)) < 2: #Only allows adjacents
                neighbours += 1 
        neighbours -= 1
        neighbours_list += str(neighbours)
    return neighbours_list
            
#----------------------------------------------------------------------------
def isNear(hits, newHit):
    
        # Finds out if a hit (newHit) is touching the cluster
        # Changing 1 to 2 would allow missing pixels for instance
        
    for hit in hits:
         if abs(hit.i-newHit.i) <= 1 and abs(hit.j-newHit.j) <= 1:
            return True

    return False

#---------------------------------------------------------------------------

def total_sum(hits):
    total = 0
    for hit in hits:
        total += hit.val
    return total    


5#--------------------------------------------------------------------------                   

def get_clusters(hits):
    
    # Makes clusters out of all adjacent hits
    
    clusters = []
    neighbours_list = get_neighbours(hits)
    i = 0
    for hit in hits:
        hit.neighbours =(neighbours_list[i])
        i += 1
                                                        #Adds hit.neighbours and hit.adjacent_neighbours to the hit
    clusters = []                   
    adjacent_neighbours_list = get_adjacent_neighbours(hits) #Could use enumerate() here?
    i = 0
    for hit in hits:                
        hit.adjacent_neighbours =(adjacent_neighbours_list[i])
        i += 1

        
    # Loop over all hits making the first that isn't already clustered
    # a new cluster
    for hit in hits:
        if hit.isClustered(): continue
        cluster = Cluster()
        cluster.add(hit)
        clusters.append(cluster)
        # Keep looping till we don't find anything
        found = True    # Necessary to start it off
        while (found):
            found = False
            # Now loop round all the other hits to see if they are touching
            for nhit in hits:
                if nhit.isClustered(): continue
                if cluster.isNear(nhit):
                    cluster.add(nhit)
                    found = True    # Makes it go around again
    return clusters
#---------------------------------------------------------------------------
class Hit:

    # Class representing a single hit

    def __init__ (self, i, j, val, neighbours, adjacent_neighbours):
    
        # Make hit from the coordinates and value
        
        self.i = i
        self.j = j
        self.val = int(val)
        self.cluster = None
        self.neighbours = neighbours
        self.adjacent_neighbours = adjacent_neighbours
    

#-----------------------------------------------------------------

    def string(self):
    
        # Make a printable string
        
        return 'i ' + str(self.i) + ' j ' + str(self.j) + ' val ' + str(self.val)

#-----------------------------------------------------------------

    def isClustered(self):
    
        # Says whether the hit is part of a cluster
        
        if self.cluster != None:
            return True
        else:
            return False

#=================================================================

class Cluster:

    # Class representing a single cluster (one or more adjacent hits)
    
    def __init__ (self):
    
        # Initialises a cluster
        
        self.hits = []
        self.sum = 0    # Sum of the values of the hits
        self.num = 0    # Number of hits
        self.x = 0      # Coordinates of centre of gravity
        self.y = 0
        self.topx = 255         # Coordinates of top edge of cluster
        self.topy = 255
        self.bottomx = 0        # Coordinates of bottom edge of cluster
        self.bottomy = 0
        self.difX = 0           # Difference in the edeges of cluster
        self.difY = 0
        self.theta = 0.0        # Angle made with horizontal by track
        self.width = 0.0
        self.height = 0.0
        self.rux = 0
        self.ruy = 255
        self.lbx = 255
        self.lby = 0
        self.difXX = 0
        self.difYY = 0
        self.roundness = 0.0
        self.hitdensity = 0.0
        self.energydensity = 0.0
        self.eph = 0.0
        self.neighbours = 0
        self.meanexcess = 0
        self.highenergyadjacent = 0
        self.type = 0
        self.energy = 0.0
#-----------------------------------------------------------------

    def add(self, hit):

        
    
        # Adds a hit to the cluster
        
        self.hits.append(hit)
        self.sum += hit.val
        self.num = len(self.hits)
        hit.cluster = self
        self.getCentre()
        self.getTopEdge()
        self.getBottomEdge()
        self.getWidth()
        self.getCorners()
        self.getHeight()
        self.getRoundness()
        self.getHitDensity()
        self.getEnergyDensity()
        self.getEPH()
        self.getNeighbours()
        self.getMeanExcess()
        self.getHighEnergyAdjacent()
        self.getType()
        self.getEquivalentEnergy()
        
#-----------------------------------------------------------------
    
    def isNear(self, newHit):
    
        # Finds out if a hit (newHit) is touching the cluster
        # Changing 1 to 2 would allow missing pixels for instance
        
        for hit in self.hits:
            if abs(hit.i-newHit.i) <= 1 and abs(hit.j-newHit.j) <= 1:
                return True

        return False
                
#-----------------------------------------------------------------

    def getCentre(self):
    
        # Crudely calculates the centre of gravity in x, y
        
        x = 0.0
        y = 0.0
        s = 0.0
        for hit in self.hits:
            x += hit.i*hit.val
            y += hit.j*hit.val
            s += hit.val
            
        self.x = int(x/s + 0.5) # Round to nearest integer
        self.y = int(y/s + 0.5)

#-----------------------------------------------------------------

    def getTopEdge(self):

                for hit in self.hits:
                        if hit.i <= self.topx and hit.j <= self.topy:
                                self.topx = hit.i
                                self.topy = hit.j

#-----------------------------------------------------------------

    def getBottomEdge(self):
            for hit in self.hits:
                        if hit.i >= self.bottomx and hit.j >= self.bottomy:
                                self.bottomx = hit.i
                                self.bottomy = hit.j
                                self.difX = self.bottomx - self.topx 
                                self.difY = self.bottomy - self.topy
                                if self.difX >= 0 and self.difY >= 0:
                                        self.theta = math.atan2(self.difY, self.difX)
                                elif self.difX == 0 and self.difY == 0:
                                        self.theta = 0.0
                                else:
                                        self.theta = (math.pi / 2.0)
                                 
#-----------------------------------------------------------------

    def getWidth(self):
                for hit in self.hits:
                        
                        if self.difX == 0 and self.difY == 0:
                                self.width = 1
                        else:
                                self.width = ((self.bottomx*math.cos(self.theta)+self.bottomy*math.sin(self.theta)) - (self.topx*math.cos(self.theta)+self.topy*math.sin(self.theta)))

#-----------------------------------------------------------------

    def getCorners(self):
            for hit in self.hits:
                        if hit.i >= self.rux and hit.j <= self.ruy:
                                self.rux = hit.i
                                self.ruy = hit.j
                        if hit.i <= self.lbx and hit.j >= self.lby:
                                self.lbx = hit.i
                                self.lby = hit.j
                                self.difXX = self.rux - self.lbx
                                self.difYY = self.lby - self.ruy

#-----------------------------------------------------------------

    def getHeight(self):

                for hit in self.hits:
                        if self.difX == 0:
                                self.height = 1.0
                        elif self.difY == 0:
                                self.height = 1.0
                        elif self.difX == 0 and self.difY == 0:
                                self.height = 1.0
                        else:
                                self.height = ((self.difXX)**2+(self.difYY)**2)**0.5

#-----------------------------------------------------------------
    def getRoundness(self):
                for hit in self.hits:
                        self.roundness = ((self.height or 1)/self.width)
                        if self.roundness > 1:
                                self.roundness = 1/self.roundness

#-----------------------------------------------------------------
    def getHitDensity(self):

                for hit in self.hits:
                        self.hitdensity = self.num/((self.width or 1)*(self.height or 1))

#-----------------------------------------------------------------
    def getEnergyDensity(self):

                for hit in self.hits:
                        self.energydensity = self.sum/((self.width or 1)*(self.height or 1))
#-----------------------------------------------------------------

    def getEPH(self):
            for hit in self.hits:
                self.eph = self.energydensity/self.hitdensity
                
         
#-----------------------------------------------------------------        

    def getNeighbours(self):
        for a in self.hits:
            if self.neighbours < a.neighbours:
                    self.neighbours = a.neighbours #Works out the highest value of Neighbours


#-----------------------------------------------------------------

    def getMeanExcess(self):
        if self.num == 1:
            self.meanexcess = 0.0
            return
        mean = float(self.sum)/float(self.num)
        selfsum = self.sum
        selfnum = float(self.num)
        for hits in self.hits:
            if float(selfsum - hits.val) / float(selfnum - 1.0) < mean: #Works out the proportion that the highest energy hit adds to the EPH
                mean = float(selfsum - hits.val) / float(selfnum - 1.0)  # i.e. it works out if the highest energy hit in the cluster is
        excessmean = (float(selfsum/selfnum) - mean)/ (selfsum/selfnum)  # much higher energy than the others, or just a bit higher than the others
        self.meanexcess = excessmean       

#-----------------------------------------------------------------

    def getHighEnergyAdjacent(self):
        i = 0
        index = 0
        value = 0
        for hit in self.hits:
            if hit.val > value:
                value = hit.val  #Works out the highest energy hit and indexes it
                index = i
            i += 1    

        self.highenergyadjacent = self.hits[index].adjacent_neighbours  #Returns the adjacent neighbours of the highest energy hit   
        



#-----------------------------------------------------------------

    def getType(self):          #Does this background flag work?
        isCoverPresent = False
        beta_coefficient = 0.0
        if self.num <= 2:       
            self.type = 'Gamma' 
            return
                
        if (isBackground and isCoverPresent): 
            if self.roundness == 1:
                beta_coefficient -= 1
            if self.neighbours == (8 or 7):
                beta_coefficient -= 1
            if self.highenergyadjacent == 4:
                beta_coefficient -= 1

            if beta_coefficient <= 2:
                self.type = 'Alpha' #Only high-energy Alphas will get through, so have to be well-defined 
                return
            else:
                self.type = 'Beta'
                return

                
        if isBackground:
            if self.num >= 10:        
                beta_coefficient += 3
                
            elif self.num >= 9:         
                beta_coefficient += 2 
                
            elif self.num >= 7:         
                beta_coefficient += 1
                                     #Note: add in width/height analysis
            if self.sum >= 200:     
                beta_coefficient += 3
                
            elif self.sum >= 125:       
                beta_coefficient += 1
                
            elif self.hitdensity >= 6:
                beta_coefficient += 1
                
            if self.roundness == 0.707106781187: #Keep this?
                beta_coefficient -= 3    

            if self.energydensity <= 10:  
                beta_coefficient += 1

            elif self.roundness <= 0.25:
                beta_coefficient += 1

            elif self.hitdensity >= 6:
                beta_coefficient += 1

            if self.highenergyadjacent == 4:
                beta_coefficient -= 3

            elif self.highenergyadjacent < 2:
                beta_coeffcient += 1

            if self.eph >= 23:
                beta_coefficient -= 1
            
            if self.meanexcess >= 0.3:
                beta_coefficient -= 2
                
            elif self.meanexcess <= 0.05:
                beta_coefficient += 1

            if self.neighbours == 8:
                beta_coefficient -= 2

            
            if beta_coefficient <= -3:
                self.type = 'Alpha'
                return                  #Obviously check these are the best ones
            elif beta_coefficient < 0:
                self.type = 'Maaybe Alpha'
                return              #Could this be going horribly wrong?
            else:
                self.type = 'Beta'
                return      #Might want to check these are being added up correctly
                



        if not(isBackground):               
            if self.neighbours == 8:
                beta_coefficient -= 2

            elif self.neighbours > 6:
                beta_coefficient -= 1
                
            if self.highenergyadjacent == 4:    
                beta_coefficient -= 1
                
            if self.roundness == 1:
                beta_coefficient -= 1

            if self.roundness == 0.707106781187:        #This section needs more work
                beta_coefficient -= 2

            if self.roundness >= 0.75:
                beta_coefficient -= 1
                
            if self.roundness < 0.2:
                beta_coefficient += 1
                
            if self.roundness < 0.1:
                beta_coefficient += 2
                
            if beta_coefficient >= 2: 
                self.type = 'Beta'     
            elif beta_coefficient <= -2:        #This section needs checking, as is simply copied/pasted
                self.type = 'Alpha'
            elif beta_coefficient == -1:
                self.type = 'Maybe Alpha'
            elif beta_coefficient == 1:
                self.type = 'Maybe Beta'
            else:
                self.type = 'Alpha or Beta'

#-----------------------------------------------------------------
    def getEquivalentEnergy(self):
        KeV = (self.sum * 0.24 + 15.78) * 1000 #Calculates KeV output 
        self.equivalentenergy = (1.6 * (10 ** -19)) * KeV   # Works out energy in joules
        if self.type == 'Alpha':
            self.equivalentenergy *= 20  #If the type is alpha, the energy value is multiplied by 20, the dose scaling factor for alphas.
#-----------------------------------------------------------------                
                
    def string(self):
    
        # Make a printable string
        
        return 'x ' + str(self.x) + ' y ' + str(self.y) + ' hits ' + str(self.num) + ' sum ' + str(self.sum)

#=================================================================

if __name__ == "__main__":
    # Program starts here

  #  outputFile = open('output.txt', "w")
    s = raw_input('Enter Directory to Process: ')

    read_in_dir(s)

    #outputFile.close()
    
    ##for file in files:
     ##   read_pix_main(file)

#=================================================================


    
