#!/usr/bin/env python3 
# Name: Your full name (CATS account username) 
# Group Members: List full names (CATS usernames) or “None”

'''
This program is designed to return calculated distances and angles between atoms in a protein crystal.

input: One line which contains coordinate information for C, N, Ca. 
output: The distance between N-C and N-Ca atoms and the angle of the C-N-Ca bond.
'''

import math
class Triad :
    """
    Calculate angles and distances among a triad of points.
 
    Author: David Bernick
    Date: March 21, 2013
    Points can be supplied in any dimensional space as long as they are consistent.
    Points are supplied as tupels in n-dimensions, and there should be three
    of those to make the triad. Each point is positionally named as p,q,r
    and the corresponding angles are then angleP, angleQ and angleR.
    Distances are given by dPQ(), dPR() and dQR()
 
    Required Modules: math
    initialized: 3 positional tuples representing Points in n-space
             p1 = Triad( p=(1,0,0), q=(0,0,0), r=(0,1,0) )
    attributes: p,q,r the 3 tuples representing points in N-space
    methods:  angleP(), angleR(), angleQ() angles measured in radians
          dPQ(), dPR(), dQR() distances in the same units of p,q,r
 
    """
 
    def __init__(self,p,q,r) :
        """ Construct a Triad. 
        
        Example construction:
            p1 = Triad( p=(1.,0.,0.), q=(0.,0.,0.), r=(0.,0.,0.) ). 
        """
        self.p = p
        self.q = q
        self.r = r
# private helper methods
    def d2 (self,a,b) : # calculate squared distance of point a to b
        return float(sum((ia-ib)*(ia-ib)  for  ia,ib in zip (a,b)))
    
    def dot (self,a,b) : # dotProd of standard vectors a,b
        return float(sum(ia*ib for ia,ib in zip(a,b)))
    
    def ndot (self,a,b,c) : # dotProd of vec. a,c standardized to b
        return float(sum((ia-ib)*(ic-ib) for ia,ib,ic in zip (a,b,c)))
    
# calculate lengths(distances) of segments PQ, PR and QR
    def dPQ (self):
        """ Provides the distance between point p and point q """
        return math.sqrt(self.d2(self.p,self.q))
    
    def dPR (self):
        """ Provides the distance between point p and point r """
        return math.sqrt(self.d2(self.p,self.r))
    
    def dQR (self):
        """ Provides the distance between point q and point r """
        return math.sqrt(self.d2(self.q,self.r))
    
    def angleP (self) :
        """ Provides the angle made at point p by segments pq and pr (radians). """
        return math.acos(self.ndot(self.q,self.p,self.r) /   math.sqrt(self.d2(self.q,self.p)*self.d2(self.r,self.p)))
    
    def angleQ (self) :
        """ Provides the angle made at point q by segments qp and qr (radians). """
        return math.acos(self.ndot(self.p,self.q,self.r) /  math.sqrt(self.d2(self.p,self.q)*self.d2(self.r,self.q)))
 
    def angleR (self) :
        """ Provides the angle made at point r by segments rp and rq (radians). """
        return math.acos(self.ndot(self.p,self.r,self.q) /  math.sqrt(self.d2(self.p,self.r)*self.d2(self.q,self.r)))

def main():
    ''' The user input will be manipulated to extract and order relevant numerical information for functional use.'''
    
    # will collect the user input as a triplet of tuples, all on one line - each tuple will have 3 elements
    C = input("Enter those tuples:")
    
    #unnecessary characters are stripped and then the tuples are split apart
    C = C.replace('Ca', "")
    C = C.replace('C', "")
    C = C.replace('N', "")
    C = C.replace('(', "")
    C = C.replace(')', "")
    C = C.split('=')
    
    #shifting over from split 
    C[0] = C[1]
    C[1] = C[2]
    C[2] = C[3]
    
    #splitting each tuple of C into 3 individual components
    C[0] = C[0].split(",")
    C[1] = C[1].split(",")
    C[2] = C[2].split(",")
    
    #type casting each component to float for later use in numerical operations
    C[0][0] = float(C[0][0])
    C[0][1] = float(C[0][1])
    C[0][2] = float(C[0][2])
    
    C[1][0] = float(C[1][0])
    C[1][1] = float(C[1][1])
    C[1][2] = float(C[1][2])
    
    C[2][0] = float(C[2][0])
    C[2][1] = float(C[2][1])
    C[2][2] = float(C[2][2])
    
    #instantiation of triad class, bearing statements assigning values to p,q,r objects 
    p1 = Triad(p = C[0], q = C[1], r = C[2])
    
    # Here are the member instantiations and print statements, answers converted to degrees and rounded. 
    distancePQ = p1.dPQ()
    print("N-C bond length = ", round(distancePQ, 2))
    
    distanceQR = p1.dQR()
    print("N-Ca bond length = ", round(distanceQR, 2))
    
    angleOfq = p1.angleQ()
    print("C-N-Ca bond angle = ", round(math.degrees(angleOfq), 1))
    
main()
