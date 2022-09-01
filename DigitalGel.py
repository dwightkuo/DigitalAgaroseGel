import sys, argparse, math, svgwrite, tempfile


### Collection of classes for drawing a Digital Gel ###


class Band:
    """Defines a Gel band


    numBases - Length of band fragment (in bp)
    annot - Annotation to display
    link - anchor link (e.g. html if the band links to something else)
    reads - the "intensity" of the band
    """
    def __init__(self, numBases, annot="", reads=1, link=None):
        self.numBases = numBases
        self.annot = annot
        self.link = link
        self.reads = reads


    def draw(self, dwg, xCoord, yCoord, intensity, dParam=None, isLadder=False):
        """Draw an individual band


        dwg - svgwrite drawing object to use
        xCoord - What x-coordinate to use to draw the band
        uCoord - What y-coordinate to use to draw the band
        intensity - brightness to draw band (can be different from # reads)
        dParam - a drawParams object that stores parameters for drawing
        """
        if dParam == None:
            dParam = drawParams()
        bandLine = dwg.rect(
                        insert = (xCoord, yCoord),
                        size = (dParam.laneWidth, 1), fill = "white", opacity = intensity )
        x = xCoord - int(float(dParam.margin)*2/3) if isLadder else 0
        y = yCoord + 3 if isLadder else dParam.gelHeight + 2 * dParam.margin - 5 


        g = dwg.g(class_="label") if not isLadder else dwg.g(class_="text")
      
        g.add( dwg.text(self.annot, insert=(x, y), 
               stroke='none', font_size='%dpx' % dParam.fontsize, 
               fill=svgwrite.rgb(100,100,100,'%')))


        if self.link:
            link = dwg.add(svgwrite.container.Hyperlink(self.link, target='_top'))
            link.add(bandLine)
            link.add(g)
            dwg.add(link)
        else:
            dwg.add(bandLine)
            dwg.add(g)    
        return


class Lane:
    """Define a Lane object
    
    name - Name of the Lane to print on top of Lane
    bands - a list of Band objects that define the lane
    isLadder - boolean to say whether this lane is a ladder
    gelHeight - the total height of the gel in pixels
    dilationFactor - the expansion factor to use when calculating band migration in the gel
    """
    def __init__(self, name, bands, dParam=None, isLadder=False):
        if dParam == None:
            dParam = drawParams()
        self.name = name
        self.bands = list(bands)
        self.isLadder = isLadder
        self.gelHeight = dParam.gelHeight
        self.dilationFactor = dParam.dilationFactor


    def __distance(self, numBases):
        """Calculate the distance a band of length 'length' would migrate in a gel given a 'dilation' factor """
        return float(self.dilationFactor) * (4 - math.log10(numBases))


    def __lengthtoIndex(self, band, low, high):
        """Convert a band distance to an index within in a list"""
        minimum = self.__distance(high)
        expFactor = self.__distance(low) - minimum


        return int((self.__distance(band.numBases) - minimum ) / 
                    float(expFactor) *(self.gelHeight-1))


    def draw(self, dwg, xCoord, minBandSize, maxBandSize, dParam=None):
        """Draw a lane
    
        dwg - svgwrite drawing object to use
        xCoord - What x-coordinate to use to draw the band
        minBandSize - minimum band size to allow
        maxBandSize - maximum band size to allow
        dParam - a drawParams object that stores parameters for drawing
        """
        minBandIntensity = 0.2
        if dParam == None:
            dParam = drawParams()
        dwg.add(dwg.text(self.name, insert = (xCoord, 10), stroke='none', 
                font_size='%dpx'%dParam.fontsize, fill=svgwrite.rgb(100,100,100,'%')))
        maxNum = max(b.reads for b in self.bands)
            
        SeenY = {}
        #Draw the bands 
        for band in self.bands:
            yCoord = self.__lengthtoIndex(band, minBandSize, maxBandSize) + dParam.margin
            intensity = 0
            if self.isLadder:
                intensity = band.reads
            else: # Normalize the band intensity (between minBandIntensity and 1)
                val = 1 - (maxNum - float(band.reads)) / float(maxNum)
                intensity = max(val, minBandIntensity)
            if yCoord not in SeenY: # never draw a band twice. Only the brightest band will be drawn
                SeenY[yCoord] = 1
                band.draw(dwg, xCoord, yCoord, intensity, dParam, self.isLadder)
        return


class Gel:
    """Gel Class
    lanes - a list of Lanes which define the Gel
    dParam - a drawParam object storing parameters for drawing
    ladderBands - where possible bands of the ladder should be drawn
    boldLadder - bands of the ladder which should be brighter
    filename - output filename for gel image
    """
    def __init__(self, lanes, dParam=None, 
        ladder = [25,50,75,100,125,150,200,250,300,400,500,600,700,800,900,1000,1250,1500,1750,2000],
#        ladder = [*range(25,150,25),*range(150,300,50),*range(300,1000,100),*range(1000,2000,250)],
        boldLadder=[50,100,200,500,1000,1500],
        filename=None):


        if dParam == None:
            dParam = drawParams()
        self.boldLadder = boldLadder
        self.lanes = list(lanes)
        self.dParam=dParam
        self.tempFile = None
        self.filename = None
        if filename:
            self.filename = filename
            self.dwg = svgwrite.Drawing(self.filename, profile='tiny')
        else:
            self.tempFile = tempfile.NamedTemporaryFile(suffix=".svg")
            self.dwg = svgwrite.Drawing(self.tempFile.name, profile='tiny')
        
        self.dwg.add_stylesheet('hover.css', title="label")
        self.maxNBases = 0
        self.minNBases = 100000
        bands = self.__makeLadder(ladder, boldLadder)
        self.ladder = Lane("Ladder", bands, dParam, isLadder=True)
    
    def __makeLadder(self, ladderBands, boldLadder):
        """Create the ladder """
        self.maxNBases = max( max(band.numBases for band in lane.bands) for lane in self.lanes)
        self.minNBases = min( min(band.numBases for band in lane.bands) for lane in self.lanes)


        bands = []
        newMin = 1000
        newMax = 0
        for ind, length in enumerate(ladderBands):
            intensity = 1.0 if length in boldLadder else 0.5
            ind1 = max(0, ind-1)
            ind2 = min(len(ladderBands)-1, ind+1)
            if (length < self.maxNBases and length > self.minNBases or
                ladderBands[ind1] < self.maxNBases and length > self.minNBases or
                ladderBands[ind2] > self.minNBases and length < self.maxNBases):
                bands.append(Band(length, annot=str(length).rjust(3,' '), reads=intensity, link=None))
                if ladderBands[ind1] < self.maxNBases and length > self.minNBases and ladderBands[ind1] < newMin:
                    newMin = ladderBands[ind1]
                if ladderBands[ind2] > self.minNBases and length < self.maxNBases and ladderBands[ind2] > newMax:
                    newMax = ladderBands[ind2]
        if newMin != 0 and newMin < self.minNBases:
            self.minNBases = newMin
        if newMax != 0 and newMax > self.maxNBases:
            self.maxNBases = newMax
        return bands


    def draw(self):
        """Draw the Gel"""
        gelWidth = self.dParam.margin * 2 + (len(self.lanes) - 1) * self.dParam.gapWidth + len(self.lanes) * self.dParam.laneWidth
        self.dwg.add(self.dwg.rect((0, 0), (gelWidth, self.dParam.gelHeight+2*self.dParam.margin), fill='black') )
        self.ladder.draw(self.dwg, self.dParam.margin, self.minNBases, self.maxNBases, self.dParam)
        for index, lane in enumerate(self.lanes):
            xCoord = self.dParam.margin + (index + 1) * (self.dParam.gapWidth + self.dParam.laneWidth)
            lane.draw(self.dwg, xCoord, self.minNBases, self.maxNBases, self.dParam)
        if self.filename:
            self.dwg.save()
        else:
            return self.dwg.tostring()
 
class drawParams:
    """ Class to store drawing parameters
    
    gelHeight - total height of gel in pixels
    margin - the edge margins of the gel (distance from any edge of the image to the gel
    laneWidth - how wide a single lane of the gel should be
    gapWidth - the distance between consecutive lanes
    fontsize - the font size for annotations
    dilationFactor - the dilation factor for gel migration
    """
    def __init__(self, gelHeight=400, margin=60, laneWidth=30, gapWidth=20, fontsize=12, dilationFactor=5):
        self.gelHeight = gelHeight 
        self.margin = margin
        self.laneWidth = laneWidth
        self.gapWidth = gapWidth
        self.fontsize =  fontsize
        self.dilationFactor = dilationFactor


if __name__ == "__main__": # Example Gel (output in svg format)


    bands2=[Band(177, reads=1925, annot="Cluster 1: 177bp | 1,925 reads", link="http://"),
        Band(51, reads=1616, annot="Cluster 2: 51bp | 1,616 reads", link="http://"),
        Band(55, reads=1327, annot="Cluster 3: 55bp | 1,327 reads", link="http://"),
        Band(154, reads=1099, annot="Cluster 4: 154bp | 1,099 reads", link="http://"),
        Band(85, reads=786, annot="Cluster 5: 85bp | 786 reads", link="http://")
        ]
    bands3=[Band(1772, reads=1925, annot="Cluster 1: 1772bp | 1,925 reads", link="http://"),
        Band(511, reads=1616, annot="Cluster 2: 551bp | 1,616 reads", link="http://"),
        Band(552, reads=1327, annot="Cluster 3: 552bp | 1,327 reads", link="http://"),
        Band(1541, reads=1099, annot="Cluster 4: 1541bp | 1,099 reads", link="http://"),
        Band(851, reads=786, annot="Cluster 5: 851bp | 786 reads", link="http://")
        ]
    Lanes = [Lane(name="1",bands = bands2), Lane(name="2",bands = bands3)]


    newGel = Gel(Lanes)
    temp = newGel.draw()
    print(temp)