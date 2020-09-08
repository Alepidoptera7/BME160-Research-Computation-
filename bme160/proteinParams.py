#!/usr/bin/env python3
# Name: Quin Lamothe (alamothe)
# Group Members: None

class ProteinParam:
    # These tables are for calculating:
    #     molecular weight (aa2mw), along with the mol. weight of H2O (mwH2O)
    #     absorbance at 280 nm (aa2abs280)
    #     pKa of positively charged Amino Acids (aa2chargePos)
    #     pKa of negatively charged Amino acids (aa2chargeNeg)
    #     and the constants aaNterm and aaCterm for pKa of the respective termini
    #  Feel free to move these to appropriate methods as you like

    # As written, these are accessed as class attributes, for example:
    # ProteinParam.aa2mw['A'] or ProteinParam.mwH2O

    aa2mw = {
        'A': 89.093, 'G': 75.067, 'M': 149.211, 'S': 105.093, 'C': 121.158,
        'H': 155.155, 'N': 132.118, 'T': 119.119, 'D': 133.103, 'I': 131.173,
        'P': 115.131, 'V': 117.146, 'E': 147.129, 'K': 146.188, 'Q': 146.145,
        'W': 204.225, 'F': 165.189, 'L': 131.173, 'R': 174.201, 'Y': 181.189
    }

    mwH2O = 18.015
    aa2abs280 = {'Y': 1490, 'W': 5500, 'C': 125}

    aa2chargePos = {'K': 10.5, 'R': 12.4, 'H': 6}
    aa2chargeNeg = {'D': 3.86, 'E': 4.25, 'C': 8.33, 'Y': 10}
    aaNterm = 9.69
    aaCterm = 2.34

    def __init__(self, inString):
        """will hold AA Composition Dictionary"""
        # replacing every keyboard symbol with ''
        inString = inString.replace('@', '')
        inString = inString.replace('~', '')
        inString = inString.replace('!', '')
        inString = inString.replace('#', '')
        inString = inString.replace('$', '')
        inString = inString.replace('%', '')
        inString = inString.replace('^', '')
        inString = inString.replace('&', '')
        inString = inString.replace('*', '')
        inString = inString.replace('(', '')
        inString = inString.replace(')', '')
        inString = inString.replace(',', '')
        inString = inString.replace('.', '')
        inString = inString.replace('<', '')
        inString = inString.replace('>', '')
        inString = inString.replace('[', '')
        inString = inString.replace(']', '')
        inString = inString.replace('=', '')
        inString = inString.replace(';', '')
        inString = inString.replace('}', '')
        inString = inString.replace('{', '')
        inString = inString.replace('|', '')
        inString = inString.replace('"', '')
        inString = inString.replace('^', '')
        inString = inString.replace('?', '')
        inString = inString.replace('1', '')
        inString = inString.replace('2', '')
        inString = inString.replace('3', '')
        inString = inString.replace('4', '')
        inString = inString.replace('5', '')
        inString = inString.replace('6', '')
        inString = inString.replace('7', '')
        inString = inString.replace('8', '')
        inString = inString.replace('9', '')
        inString = inString.replace('0', '')
        inString = inString.upper()

        self.initAAcomp = {
            'A': inString.count('A'), 'G': inString.count('G'),
            'M': inString.count('M'), 'S': inString.count('S'),
            'C': inString.count('C'), 'H': inString.count('H'),
            'N': inString.count('N'), 'T': inString.count('T'),
            'D': inString.count('D'), 'I': inString.count('I'),
            'P': inString.count('P'), 'V': inString.count('V'),
            'E': inString.count('E'), 'K': inString.count('K'),
            'Q': inString.count('Q'), 'W': inString.count('W'),
            'F': inString.count('F'), 'L': inString.count('L'),
            'R': inString.count('R'), 'Y': inString.count('Y')
        }

        # multipliers for + charges
        self.kCount = self.initAAcomp['K']
        self.rCount = self.initAAcomp['R']
        self.hCount = self.initAAcomp['H']

        # multipliers for - charges
        self.dCount = self.initAAcomp['D']
        self.eCount = self.initAAcomp['E']
        self.cCount = self.initAAcomp['C']
        self.yCount = self.initAAcomp['Y']

        # + charge pka values
        kCharge = self.aa2chargePos['K']
        rCharge = self.aa2chargePos['R']
        hCharge = self.aa2chargePos['H']

        # - charge pka values
        dCharge = self.aa2chargeNeg['D']
        eCharge = self.aa2chargeNeg['E']
        cCharge = self.aa2chargeNeg['C']
        yCharge = self.aa2chargeNeg['Y']

        # for use in excluding the charges of aminos from the pI calculation
        self.pospKaVal = []

        self.pospKaVal.insert(0, kCharge)
        self.pospKaVal.insert(1, rCharge)
        self.pospKaVal.insert(2, hCharge)

        # - charge presence loads list with the peptides - charged amino values
        self.negpKaVal = []

        self.negpKaVal.insert(0, dCharge)
        self.negpKaVal.insert(1, eCharge)
        self.negpKaVal.insert(2, cCharge)
        self.negpKaVal.insert(3, yCharge)


    def aaCount(self):
        """Will count the number of occurrences of each amino acid in the
        sequence and form a sum of the values."""

        self.sumAA = sum(self.initAAcomp.values())

        self.sumAA = int(self.sumAA)

        return int(self.sumAA)

    def aaComposition(self):
        """Will return the AA composition dictionary of init method as a list."""

        newHolder = self.initAAcomp

        return newHolder

    def _charge_(self, pH):
        """Zill aid in finding the pH of the given peptide's Zwitterion."""
        # forming + charge quotient
        # each positive amino charge value is evaluated separately in an attempt to reduce error, still returning the same number!

        negTotal = 0
        posTotal = 0

        #Dennis provided an alteration to my method which caused it to look like this! Better use of dictionaries.

        for aminoAcid, pKa in self.aa2chargePos.items():
            posQuotientNumerator = (10 ** pKa)
            posQuotientDenominator = (10 ** pKa) * (10 ** pH)
            posTotal += self.initAAcomp[aminoAcid] * (posQuotientNumerator / posQuotientDenominator)

        posQuotientNumerator = (10 ** self.aaNterm)
        posQuotientDenominator = (10 ** self.aaNterm) * (10 ** pH)
        posTotal += (posQuotientNumerator / posQuotientDenominator)

        #and now to calculate quotients for each of the present negative charged aminos at given pH

        for aminoAcid, pKa in self.aa2chargeNeg.items():
            negQuotientNumerator = (10 ** pH)
            negQuotientDenominator = (10 ** pH) * (10 ** pKa)
            negTotal += self.initAAcomp[aminoAcid] * (negQuotientNumerator / negQuotientDenominator)

        negQuotientNumerator = (10 ** pH)
        negQuotientDenominator = (10 ** pH) * (10 ** self.aaCterm)
        negTotal += (negQuotientNumerator / negQuotientDenominator)

        netCharge = posTotal - negTotal

        return netCharge

    def pI(self):
        """The pI of the peptide will be estimated by testing the results of bisection of values between
         0.0 pH and 14 pH. The results of the charge function will be evaluated to
        determine the pH value that causes the charge of the peptide to be 0.0"""


        #Mohammad helped me with this function in office hours by use of pseudo code-- I dont have it working still.
        Low = 0.0
        High = 14

        while (High - Low) > 0.01:

            mid = (High - Low) / 2
            # If 0.0 is found at mid, this is the zwitterion
            thisCharge = self._charge_(mid)
            # If the charge is less than 0.0, the peptide is still acidic and the high pH should be lowered.
            if thisCharge < 0.0:
                mid = High

            # If the charge is more than 0.0, the peptide is basic and the low should be raised.
            else:
                Low = mid

            # If we reach here, then the element was not present
        return mid


    def molarExtinction(self):
        """Will multiply the number of Y, W, & C by the respective extinction coefficients"""

        yCount = self.initAAcomp['Y']
        wCount = self.initAAcomp['W']
        cCount = self.initAAcomp['C']

        yEx = self.aa2abs280['Y']
        wEx = self.aa2abs280['W']
        cEx = self.aa2abs280['C']

        yCountyEx = yCount * yEx
        wCountwEx = wCount * wEx
        cCountcEx = cCount * cEx

        extinction = yCountyEx + wCountwEx + cCountcEx

        return extinction

    def massExtinction(self):
        myMW = self.molecularWeight()
        return self.molarExtinction() / myMW if myMW else 0.0

    def molecularWeight(self):

        """multiply the AA occurance dictionary values by the corresponding molecular weights"""

        quantMWmulti = {k: self.initAAcomp[k] * self.aa2mw[k] for k in self.initAAcomp}

        sumQuantMWmulti = sum(quantMWmulti.values())

        sumAA = sum(self.initAAcomp.values())

        sumAA = int(sumAA)

        # the number of H2O molecules lost to dehydration would be half of the number of aminos bonded
        numberOfH2O = (sumAA - 1)

        # the number of AA -2, multiplied by mw of H20 to account for the total mw of H20
        massH2O = numberOfH2O * self.mwH2O

        sumQuantMWmulti = sumQuantMWmulti - massH2O

        return sumQuantMWmulti

    # Please do not modify any of the following.  This will produce a standard output that can be parsed


import sys


def main():
    """Collect input, remove non alphabet character, convert to upper case, repeat."""

    inString = input('protein sequence?')

    while inString:
        myParamMaker = ProteinParam(inString)
        myAAnumber = myParamMaker.aaCount()
        print("Number of Amino Acids: {aaNum}".format(aaNum=myAAnumber))
        print("Molecular Weight: {:.1f}".format(myParamMaker.molecularWeight()))
        print("molar Extinction coefficient: {:.2f}".format(myParamMaker.molarExtinction()))
        print("mass Extinction coefficient: {:.2f}".format(myParamMaker.massExtinction()))
        print("Theoretical pI: {:.2f}".format(myParamMaker.pI()))
        print("Amino acid composition:")

        myAAcomposition = myParamMaker.aaComposition()
        keys = list(myAAcomposition.keys())
        keys.sort()

        if myAAnumber == 0: myAAnumber = 1  # handles the case where no AA are present
        for key in keys:
            print("\t{} = {:.2%}".format(key, myAAcomposition[key] / myAAnumber))

        inString = input('protein sequence?')

        # VLSPADKTNVKAAW


if __name__ == "__main__":
    main()
