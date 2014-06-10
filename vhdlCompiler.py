'''
CS232 Project 5 extension
Due on Mar 21, 2014
@author: Ian Tibbetts and Ryan Newell
'''
import sys
#found this implementation at 
#http://stackoverflow.com/questions/36932/how-can-i-represent-an-enum-in-python
def enum(*sequential, **named):
    enums = dict(zip(sequential, range(len(sequential))), **named)
    return type('Enum', (), enums)

class vhdlCompiler:
    def __init__(self):
        self.state = enum("START",
                          "LOOP",
                          "ERROR")
        
        self.keywords = ["ON","OFF","SHIFT_LEFT","SHIFT_RIGHT",
                         "ROTATE_LEFT","ROTATE_RIGHT","INVERT"]
        self.keywordDefinitions = ("keywords:\n" +
                                   "Do_# - begin a loop that iterates # times\n" +
                                   "Loop - end a loop\n" +
                                   "Set_######## - sets each of the 8 lights to their #, either 1 or 0.\n" +
                                   "On - turn all 8 lights on\n" +
                                   "Off - turn all 8 lights off\n" +
                                   "Shift_Left - shift 8 lights left by one\n" +
                                   "Shift_Right - shift 8 lights right by one\n" +
                                   "Rotate_Left - rotate rightmost of 8 lights to the leftmost\n" +
                                   "Rotate_Right - rotate leftmost of 8 lights to the rightmost\n" +
                                   "Invert - flip all on lights to off and visa versa")
        self.setBits = ""
        self.errorMessage = ""
        self.startLoopAddr = 0
        self.address = 0
        self.addrLength = 8
        self.instrLength = 10
        self.loopCount = 0
    
    def parseSource(self,sourceText):
        outputStr = "data <= \n"
        curState = self.state.START
        for word in sourceText.split():
            word = word.upper()
            stateStr = self.getStateStr(curState)
            print("In state: {}".format(stateStr),end=" ")
            nextState = self.getNextState(word, curState)
            stateStr = self.getStateStr(nextState)
            print("read word: \"" + word +
                  "\" moved to state: {}".format(stateStr))
            if nextState == self.state.ERROR:
                return ""
            else:
                outputStr += self.getMachineInstructions(word)
            curState = nextState
            
        outputStr += self.getMachineInstructions("")
        return outputStr
          
    def getMachineInstructions(self,word):
        instrList = []
        instructions = "-- " + word + "\n"
        if word.split("_")[0] == "DO":
            loopCountStr = "{0:b}".format(self.loopCount)
            if len(loopCountStr) > self.addrLength:
                self.errorMessage += ("Error: new address exceeds give address length" +
                                      " of {}\n".format(self.addrLength))
            while len(loopCountStr) < self.addrLength:
                loopCountStr = "0" + loopCountStr
            instrList += [["\"001010" + loopCountStr[4:8] + "\"","set low 4 bits of ACC"]]
            instrList += [["\"001110" + loopCountStr[0:4] + "\"","set high 4 bits of ACC"]]
            instrList += [["\"0110000100\"","move ACC into LOOP"]]
        elif word == "LOOP":
            instructions += self.subOneFromLOOP()
            branchZWhenAddr = self.getNewAddress()
            branchUWhenAddr = self.getNewAddress()
            branchZToAddr = self.getNewAddress()
            instructions += ("\"11" + branchZToAddr + "\"" + 
                             " when addr = \"" + 
                             branchZWhenAddr +
                             "\" else -- " +
                             "break when LOOP is zero\n")
            instructions += ("\"10" + self.getStartLoopAddr() + "\""
                             " when addr = \"" + 
                             branchUWhenAddr +
                             "\" else -- " +
                             "branch unconditional to top of loop otherwise\n")
        elif word == "ON":
            instrList += [["\"0001110000\"","move all 1's into LR"]]
        elif word == "OFF":
            instrList += [["\"0000110000\"","move all 1's to ACC"]]
            instrList += [["\"0110011000\"","xor ACC with all 1's"]]
            instrList += [["\"0001000000\"","move ACC to LR"]]
        elif word == "SHIFT_LEFT":
            instrList += [["\"0000010000\"","move LR to ACC"]]
            instrList += [["\"0101000000\"","shift ACC left"]]
            instrList += [["\"0001000000\"","move ACC to LR"]]
        elif word == "SHIF[T_RIGHT":
            instrList += [["\"0000010000\"","move LR to ACC"]]
            instrList += [["\"0101100000\"","shift ACC right"]]
            instrList += [["\"0001000000\"","move ACC to LR"]]
        elif word == "ROTATE_LEFT":
            instrList += [["\"0000010000\"","move LR to ACC"]]
            instrList += [["\"0111000000\"","rotate ACC left"]]
            instrList += [["\"0001000000\"","move ACC to LR"]]
        elif word == "ROTATE_RIGHT":
            instrList += [["\"0000010000\"","move LR to ACC"]]
            instrList += [["\"0111100000\"","rotate ACC right"]]
            instrList += [["\"0001000000\"","move ACC to LR"]]
        elif word == "INVERT":
            instrList += [["\"0000010000\"","move LR to ACC"]]
            instrList += [["\"0110011000\"","xor ACC with all 1's"]]
            instrList += [["\"0001000000\"","move ACC to LR"]]
        elif word.split("_")[0] == "SET":
            instrList += [["\"001010" + self.setBits[4:8] + "\"","set low 4 bits of ACC"]]
            instrList += [["\"001110" + self.setBits[0:4] + "\"","set high 4 bits of ACC"]]
            instrList += [["\"0001000000\"","move ACC to LR"]]
        else:
            return ("\"1000000000\"; -- branch to beginning when no more instructions\n")
        for instr in instrList:
            newAddr = self.getNewAddress()
            instructions += (instr[0] + 
                             " when addr = \"" + 
                             newAddr +
                             "\" else -- " +
                             instr[1] + "\n")
        if word.split("_")[0] == "DO":
            self.startLoopAddr = self.address
        return instructions
    
    def getNextState(self,word, curState):
        if curState == self.state.START:
            if word.split("_")[0] == "DO":
                if word.split("_")[-1].isdigit():
                    self.loopCount = int(word.split("_")[-1])
                    curState = self.state.LOOP
                else:
                    self.errorMessage += ("Error: DO_# - Do must have number after an _\n")
                    curState = self.state.ERROR
            elif word == "LOOP":
                self.errorMessage += ("Error: \"Loop\" without \"Do\"\n")
                curState = self.state.ERROR
            elif word.split("_")[0] == "SET":
                if word.split("_")[-1].isdigit():
                    bitCount = 0
                    for bit in word.split("_")[-1]:
                        if not(bit == "0" or bit == "1"):
                            self.errorMessage += ("Error: SET_######## - bit {} is invalid\n".format(bit))
                        bitCount += 1
                    if bitCount != 8:
                        self.errorMessage += ("Error: SET_######## - Do must 8 0's or 1's after an _\n")
                    self.setBits = word.split("_")[-1]
                    curState = self.state.START
                else:
                    self.errorMessage += ("Error: SET_######## - Do must 8 0's or 1's after an _\n")
                    curState = self.state.ERROR
            elif self.keywords.count(word):
                curState = self.state.START
            else:
                self.errorMessage += ("Error: word \"{}\" not recognized" + 
                                      " only use {}\n"
                                      ).format(word,
                                               self.keywordDefinitions)
                curState = self.state.ERROR
                
        elif curState == self.state.LOOP:
            if word.split("_")[0] == "DO":
                self.errorMessage = ("Error: \"Do\" inside loop.\n" +
                                     " do note support nested loops.\n")
                curState = self.state.ERROR
            elif word == "LOOP":
                curState = self.state.START
            elif word.split("_")[0] == "SET":
                if word.split("_")[-1].isdigit():
                    bitCount = 0
                    for bit in word.split("_")[-1]:
                        if not(bit == "0" or bit == "1"):
                            self.errorMessage += ("Error: SET_######## - bit {} is invalid\n".format(bit))
                        bitCount += 1
                    if bitCount != 8:
                        self.errorMessage += ("Error: SET_######## - Do must 8 0's or 1's after an _\n")
                    self.setBits = word.split("_")[-1]
                    curState = self.state.LOOP
                else:
                    self.errorMessage += ("Error: SET_######## - Do must 8 0's or 1's after an _\n")
                    curState = self.state.ERROR
            elif self.keywords.count(word):
                curState = self.state.LOOP
            else:
                self.errorMessage += ("Error: word \"{}\" not recognized" + 
                                      " only use {}\n"
                                      ).format(word,
                                              self.keywordDefinitions)
                curState = self.state.ERROR
        else: #state.ERROR
            curState = self.state.ERROR
            
        return curState
                
    def getStateStr(self,curState):
        if curState == self.state.START:
            return "START"
        elif curState == self.state.LOOP:
            return "LOOP"
        else: #ERROR
            return "ERROR"
        
    def addOneToLOOP(self):
        return ("\"0100010101\"" + 
                " when addr = \"" + 
                self.getNewAddress() +
                "\" else -- " +
                "increment loop register by 1\n")
            
    def subOneFromLOOP(self):
        return ("\"0100110101\"" + 
                " when addr = \"" + 
                self.getNewAddress() +
                "\" else -- " +
                "decrement loop register by 1\n")
        
    def getNewAddress(self):
        newAddr = "{0:b}".format(self.address)
        if len(newAddr) > self.addrLength:
            self.errorMessage += ("Error: new address exceeds give address length" +
                                  " of {}\n".format(self.addrLength))
        while len(newAddr) < self.addrLength:
            newAddr = "0" + newAddr
        self.address += 1
        return newAddr
    
    def getStartLoopAddr(self):
        startAddr = "{0:b}".format(self.startLoopAddr)
        if len(startAddr) > self.addrLength:
            self.errorMessage += ("Error: start address exceeds give address length" +
                                  " of {}\n".format(self.addrLength))
        while len(startAddr) < self.addrLength:
            startAddr = "0" + startAddr
        self.address += 1
        return startAddr
        
if __name__ == "__main__":
    print(sys.argv)
    print(len(sys.argv))
    compiler = vhdlCompiler()
    if len(sys.argv) < 2 or len(sys.argv) > 3:
        print("usage: expect input file name as command line input.\n" +
              "Output file is an optional second parameter.\n" +
              "Input file should be program instruction written in" +
              " plain english using {}".format(compiler.keywordDefinitions))
    else:
        print("input file is " + sys.argv[1])
        sourceFile = open(sys.argv[1], encoding ="utf-8")
        sourceText = sourceFile.read()
        sourceFile.close()
        print("Input File Contents:\n" + sourceText)
        print(sourceText.split())
        
        parseResult = compiler.parseSource(sourceText)
        if compiler.errorMessage:#if there is an error string
            print(parseResult)
            print("Error(s) occurred:\n{}".format(compiler.errorMessage))
        else:
            print(parseResult)
            if len(sys.argv) == 2:
                listOutputPath = sys.argv[0].split("\\")[:-1]
                listOutputFile = listOutputPath + ["outputROM.vhd"]
                outputFileName = "outputROM"
                print("output file is " + 
                      "\\".join(listOutputFile))
                outputFile = open("\\".join(listOutputFile), mode="w", encoding ="utf-8")
            elif len(sys.argv) == 3:
                print("output file is " + sys.argv[2])
                outputFileName = sys.argv[2].split("\\")[-1].split(".")[0]
                outputFile = open(sys.argv[2], mode="w", encoding ="utf-8")
            else:
                print("usage: expect input file name as command line input.\n" +
                      "Output file is an optional second parameter.\n" +
                      "Input file should be program instruction written in" +
                      " plain english using {}".format(compiler.keywordDefinitions))
                      
            outputStr = "-- Ian Tibbetts and Ryan Newell\n"
            outputStr += "-- Project 5 " + outputFileName + ".vhd\n"
            outputStr += "-- Due: Mar 21, 2014\n"
            outputStr += "library ieee;\n"
            outputStr += "use ieee.std_logic_1164.all;\n"
            outputStr += "use ieee.numeric_std.all;\n"
            outputStr += "\n"
            outputStr += "entity " + outputFileName + " is\n"
            outputStr += "  port (\n"
            outputStr += "    addr : in std_logic_vector (7 downto 0);\n"
            outputStr += "    data : out std_logic_vector (9 downto 0));\n"
            outputStr += "end entity;\n"
            outputStr += "\n"
            outputStr += "architecture rtl of " + outputFileName + " is\n"
            outputStr += "\n"
            outputStr += "begin\n"
            outputStr += "\n"
            outputStr += parseResult
            outputStr += "\n"
            outputStr += "end rtl;\n"
            outputFile.write(outputStr)
            outputFile.close()