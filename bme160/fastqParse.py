#!/usr/bin/env python3 
# Name: Quin Lamothe (alamothe)
# Group Members: None 

'''
This program is designed to parse a FASTQ label to extract and display useful information to the user. 

input: FASTQ label in one line. 
output: Seven individual elements of information. 
'''

class FastqString (str):
    ''' The user input will be parsed and relevant information will be extracted and displayed.'''
    def parse(self):
        ''' This is the member that will strip characters, split the string, and display the information.'''
        #user input collection
        fastQ = input("enter FASTQ sequence.")
        
        #Removing unnecessary characters and splitting the string into tokens. 
        fastQ = fastQ.replace("@","")
        fastQ = fastQ.replace("_","")
        fastQ = fastQ.replace("*","")
        fastQ = fastQ.split(":")
        #Printing each token with corresponding context statement. 
        print("This is it:")
        print( "instrument = ", fastQ[0])
        print("Run ID = ", fastQ[1])
        print("Flow Cell ID = ", fastQ[2])
        print("Flow Cell Lane = ", fastQ[3])
        print("Tile Number = ", fastQ[4])
        print("X-coord = ", fastQ[5])
        print("Y-coord = ", fastQ[6])
        
def main():
    ''' Istantiation and call of the class which will do all of the necessary work to extract and display information.'''

main()
#class instantiation and call to member function 
fastQ = str 
thisFast = FastqString(fastQ)
outPut = thisFast.parse()