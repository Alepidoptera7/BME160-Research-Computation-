#!/usr/bin/env python3 
# Name: Quin Lamothe (alamothe)
# Group Members: None 

'''
Read a DNA string from user input and return a collapsed substring of embedded Ns to: {count}.

Example: 
input: AaNNNNNNGTC
output: AA{6}GTC

Any lower case letters are converted to uppercase
'''

class DNAstring (str):
    def length (self):
        return (length(self))
    
    def purify(self):
        ''' Return an upcased version of the string, collapsing a single run of Ns.'''
        
        # this is my host of variable definitions 
        
        data = input("some DNA?")
        data = data.upper()
        
        #the position of the first 'N' is found in the string, also the number of 'N's
        firstN = data.find("N")
        countN = data.count("N")
        
        #then the string is split starting at the position in the string where the first N was found
        dataNew = data.split(data[firstN], countN)
        #type casting to string again, after splitting 
        countN = str(countN)
        #manufacturing that set of curly brackets that will hold the number of 'N's in the final output 
        curlyThings = "{" +  countN + "}"      
        
        #the final output is assembled 
        data = dataNew[0] + curlyThings + dataNew[-1]
       
        return(data)
    
def main():
    ''' Get user DNA data and clean it up.'''
    
    data = str
    #class instantiations and function call 
    thisDNA = DNAstring (data)
    pureData = thisDNA.purify()
    print (pureData)
    
main()