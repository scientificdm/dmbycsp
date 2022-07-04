import numpy as np
import sys
import itertools
import timeit

fileName = sys.argv[1] # 'chembl'
fileNameTr = sys.argv[2] # 'ABL_CHEMBL_10percent_support30.csv'
fileNameAct = sys.argv[3] # 'ABL1_ChEMBL_NV_MW_800_percent10.txt'

THETA = int(sys.argv[4]) # 32 -- chi-square threshold
THETA2 = int(sys.argv[5]) # 8 -- max pattern size
THETA3 = int(sys.argv[6]) # 20 -- minimum frequency

class Variable():
    def __init__(self,label,relativeIndex,rootIndex):
        self.name = label
        self.type = label[0]
        self.index = rootIndex
        self.relative_index = relativeIndex
    
    def getName(self):
        return self.name
    
    def getType(self):
        return self.type
    
    def getIndex(self):
        return self.index
    
    def getRelativeIndex(self):
        return self.relative_index
        
##----------------------------------
class Domain():
    def __init__(self,domainValues):
        self.values = domainValues
        self.eliminated = []
        
    def getValues(self):
        return self.values
    
    def getMin(self):
        return min(self.values)

    def getMax(self):
        return max(self.values)
    
    def getLen(self):
        return len(self.values)
    
    def isEmpty(self):
        return len(self.values) == 0
        #return self.values == []
        
    def removeValue(self,value):
        self.values.remove(value)
        self.eliminated.append(value)
        
    def addValue(self,value):
        self.values.append(value)
        self.eliminated.remove(value)

##----------------------------------
class Constraint():
    def __init__(self,constraintId,constraintType,firstVariable,secondVariable=None,operationType=None,unaryValue=None):
        self.name = constraintId
        self.type = constraintType
        self.var1 = firstVariable
        self.var2 = secondVariable
        self.var_relative_index = firstVariable.getRelativeIndex()
        self.var_index = firstVariable.getIndex()
        self.var_type = firstVariable.getType()
        self.operation = operationType
        self.value = unaryValue

    def getName(self):
        return self.name
        
    def getType(self):
        return self.type

    def getFirstVariable(self):
        return self.var1
        
    def getSecondVariable(self):
        return self.var2
    
    def getVariableIndex(self):
        return self.var_index
    
    def getVariableRelativeIndex(self):
        return self.var_relative_index
    
    def getConstraint(self):
        return (self.name,self.type,self.var1,self.var2,self.operation,self.value)
    
    def getVariableType(self):
        return self.var_type
    
    def xsqaure(self,p,n):
        if (p+n) != 0:
            return (p-(p+n)*DPD)**2/((p+n)*DPD) + (n-(p+n)*DND)**2/((p+n)*DND) + (DP-p-(D-p-n)*DPD)**2/((D-p-n)*DPD) + (DN-n-(D-p-n)*DND)**2/((D-p-n)*DND)
        else:
            return 0
    
    def verify(self,firstValue,secondValue,orientation,items_min=None,items_max=None,transactions_min=None,transactions_max=None):
        if self.type == 'unary':
            # unary constraints:
            if (self.operation == '>') and (firstValue > self.value):
                return True
            elif (self.operation == '<') and (firstValue < self.value):
                return True
            elif (self.operation == '=') and (firstValue == self.value):
                return True
            elif (self.operation == '!=') and (firstValue != self.value):
                return True
            elif (self.operation == '<=') and (firstValue <= self.value):
                return True
            elif (self.operation == '>=') and (firstValue >= self.value):
                return True
            else:
                return False
        elif self.type == 'binary':
            # binary constraints:
            if (self.operation == '>') and (firstValue > secondValue):
                return True
            elif (self.operation == '<') and (firstValue < secondValue):
                return True
            elif (self.operation == '=') and (firstValue == secondValue):
                return True
            elif (self.operation == '!=') and (firstValue != secondValue):
                return True
            elif (self.operation == '<=') and (firstValue <= secondValue):
                return True
            elif (self.operation == '>=') and (firstValue >= secondValue):
                return True
            elif (self.operation == '>>') and (firstValue + self.value > secondValue):
                return True
            elif (self.operation == '<<') and (firstValue + self.value < secondValue):
                return True
            elif (self.operation == '==') and (firstValue + self.value == secondValue):
                return True
            elif (self.operation == '!==') and (firstValue + self.value != secondValue):
                return True
            elif (self.operation == '<<=') and (firstValue + self.value <= secondValue):
                return True
            elif (self.operation == '>>=') and (firstValue + self.value >= secondValue):
                return True
            else:
                return False
        elif self.type == 'coverage':
            # coverage constraint:
            if ((firstValue == 1) and not (np.dot(DB2[self.var_relative_index,:],items_min) > 0)) or ((firstValue == 0) and not (np.dot(DB2[self.var_relative_index,:],items_max) == 0)):
                return True
            else:
                return False
        elif self.type == 'frequency':
            # frequency constraint:
            if ((firstValue == 1) and (np.dot(DB[:,self.var_relative_index],transactions_max) >= THETA3)) or (firstValue == 0):
                  return True
            else:
                return False
        elif self.type == 'emerging':
            # emerging constraint:
            if ((firstValue == 1) and (max(self.xsqaure(np.dot(DB1[:,self.var_relative_index],transactions_max), np.dot(DB0[:,self.var_relative_index],transactions_min)), self.xsqaure(np.dot(DB1[:,self.var_relative_index],transactions_min), np.dot(DB0[:,self.var_relative_index],transactions_max))) >= THETA)) or (firstValue == 0):
                  return True
            else:
                return False
        elif self.type == 'closure':
            # closure constraint:
            if ((firstValue == 1) and not (np.dot(DB2[:,self.var_relative_index],transactions_min) > 0)) or ((firstValue == 0) and not (np.dot(DB2[:,self.var_relative_index],transactions_max) == 0)):
                  return True
            else:
                return False
        elif self.type == 'size':
            # max size constraint:
            if ((firstValue == 1) and not (np.dot(DB[self.var_relative_index,:],items_min) > THETA2)) or (firstValue == 0):
                return True
            else:
                return False
            
##----------------------------------
class CSP():
    def __init__(self,numItems,numTrans,setVariables,setDomains,setConstraints):
        self.num_items = numItems
        self.num_transactions = numTrans
        self.variables = setVariables
        self.domains = setDomains
        self.items_min = []
        self.items_max = []
        self.transactions_min = []
        self.transactions_max = []
        for d in setDomains:
            if len(self.items_min) < numItems:
                self.items_min.append(d.getMin())
                self.items_max.append(d.getMax())
            else:
                self.transactions_min.append(d.getMin())
                self.transactions_max.append(d.getMax())
        self.constraints = setConstraints
        self.max_constraint_id = 2*numTrans+3*numItems
        self.removed_step = {}
        self.not_instantiated = []
        self.num_pos_patterns = 0
        self.num_neg_patterns = 0
    
    def getNumberVariables(self):
        return len(self.variables)
    
    def getVariableName(self,domain):
        return self.variables[self.domains.index(domain)].getName()
           
    def getDomain(self,variableIndex):
        return self.domains[variableIndex]
        
    def getNumberConstraints(self):
        return len(self.constraints)
    
    def getMaxConstraintId(self):
        self.max_constraint_id+=1
        return self.max_constraint_id
    
    def chisquare(self,p,n):
        if ((p+n) != 0) and ((p+n) != D):
            return (p-(p+n)*DPD)**2/((p+n)*DPD) + (n-(p+n)*DND)**2/((p+n)*DND) + (DP-p-(D-p-n)*DPD)**2/((D-p-n)*DPD) + (DN-n-(D-p-n)*DND)**2/((D-p-n)*DND)
        else:
            return 0
        
    def removeMinMax(self,variableType,variableId,value):
        if variableType == 'I':
            if value == 0:
                if self.items_max[variableId] == 1:
                    self.items_min[variableId] = 1
                else:
                    self.items_min[variableId] = None
                    self.items_max[variableId] = None
            else:
                if self.items_min[variableId] == 0:
                    self.items_max[variableId] = 0
                else:
                    self.items_min[variableId] = None
                    self.items_max[variableId] = None
        elif variableType == 'T':
            if value == 0:
                if self.transactions_max[variableId] == 1:
                    self.transactions_min[variableId] = 1
                else:
                    self.transactions_min[variableId] = None
                    self.transactions_max[variableId] = None
            else:
                if self.transactions_min[variableId] == 0:
                    self.transactions_max[variableId] = 0
                else:
                    self.transactions_min[variableId] = None
                    self.transactions_max[variableId] = None
                    
    def restoreMinMax(self,variableType,variableId,value):
        if variableType == 'I':
            if value == 0:
                if self.items_max[variableId] == None:
                    self.items_min[variableId] = 0
                    self.items_max[variableId] = 0
                else:
                    self.items_min[variableId] = 0
            else:
                if self.items_min[variableId] == None:
                    self.items_min[variableId] = 1
                    self.items_max[variableId] = 1
                else:
                    self.items_max[variableId] = 1
        elif variableType == 'T':
            if value == 0:
                if self.transactions_max[variableId] == None:
                    self.transactions_min[variableId] = 0
                    self.transactions_max[variableId] = 0
                else:
                    self.transactions_min[variableId] = 0
            else:
                if self.transactions_min[variableId] == None:
                    self.transactions_min[variableId] = 1
                    self.transactions_max[variableId] = 1
                else:
                    self.transactions_max[variableId] = 1
            
    def removeDomainValue(self,variable,value):
        self.removeMinMax(variable.getType(),variable.getRelativeIndex(),value)
        self.getDomain(variable.getIndex()).removeValue(value)
        
    def restoreDomainValue(self,variable,value):
        self.restoreMinMax(variable.getType(),variable.getRelativeIndex(),value)
        self.getDomain(variable.getIndex()).addValue(value)
        
    def printCSP(self):
        variables = []
        for v in self.variables:
            variables.append(v.getName())
        print('Variables: '+str(variables))
        domains = []
        for d in self.domains:
            domains.append(d.getValues())
        print('Domains: '+str(domains))
        constraints = []
        for c in self.constraints:
            constraints.append(c.getConstraint())
        print('Constraints: '+str(constraints))
        
    def printDomains(self):
        domains = []
        for d in self.domains:
            domains.append(d.getValues())
        print('Domains: '+str(domains))
        
    def checkSolution(self):      
        # Determine associated positive transactions:
        num_trans_pos = np.dot(DC,self.transactions_max)
        
        # Determine associated negative transactions:
        num_trans_neg = np.dot(CD,self.transactions_max)
        
        # Verify discriminating measure:
        f_measure = self.chisquare(num_trans_pos,num_trans_neg)
        
        if f_measure >= THETA:
            # Determine class of pattern:
            if num_trans_pos > num_trans_neg:
                self.num_pos_patterns+=1
            else:
                self.num_neg_patterns+=1
                
            # Determine associated positive transactions:
            trans_pos = [str(j)+' ' for j in (self.transactions_max*DC).nonzero()[0]]
            
            # Determine associated negative transactions:
            trans_neg = [str(j)+' ' for j in (self.transactions_max*CD).nonzero()[0]]

            # Determine frequent itemset:
            itemset = [PHARMA[i]+' ' for i in range(numItems) if np.max(self.domains[i].getValues()) == 1]
            
            # Print the result:
            print(''.join(itemset).strip()+' ('+str(num_trans_pos+num_trans_neg)+': +'+str(num_trans_pos)+' -'+str(num_trans_neg)+') '+str(np.round(f_measure,2)))
            
            # Output the result to a file:
            f_out.write(str(numSolutions+1)+',"'+''.join(itemset).strip()+'","'+''.join(trans_pos).strip()+'","'+''.join(trans_neg).strip()+'",'+str(np.round(f_measure,2))+'\n')
            
            return True
        else:
            return False
        
    def mac(self):
        self.not_instantiated = []
        for i in range(self.num_items):
            self.removed_step[i+1] = []
            self.not_instantiated.append(self.variables[i])
        return (self.establishAC() and self.searchSolutionMAC(1))
    
    def establishAC(self):

        return self.initialization()
       
    def initialization(self):     
        # process unary constraints:
        for c in self.constraints:
            if c.getType() == 'unary':
                if self.reviseUnary(self.getDomain(c.getVariableIndex()),c,0):
                    if self.getDomain(c.getVariableIndex()).isEmpty():

                         return False
        
        Q = [] # queue of binary constraints to process
        for c in self.constraints:
            if c.getType() != 'unary':
                Q.append(c)
        
        while Q != []:
            c = Q.pop(0)
                        
            # process coverage constraints:
            if (c.getType() == 'coverage') or (c.getType() == 'size'):
                #!print('process coverage constraint c'+str(c.getName()))
                if self.reviseUnary(self.getDomain(c.getVariableIndex()),c,0):
                    #if self.getDomain(c.getVariableIndex()).getValues() == []:
                    if self.getDomain(c.getVariableIndex()).isEmpty():
                        # print explanation:
                        #!self.failureExplanation(c.getFirstVariable())
                    
                        return False
                    
                    #addedFlag = False
                    if self.getDomain(c.getVariableIndex()).getValues()[0] != 1:
                        # add frequency constraints:
                        freq_const_id = [i for i in range(self.num_items) if DB[c.getVariableRelativeIndex(),i] != 0]
                        for i1 in freq_const_id:
                            c1 = self.constraints[i1+self.num_transactions]
                            if 1 in self.getDomain(i1).getValues():
                                if c1 not in Q:
                                    Q.append(c1)
                                    #if not addedFlag:
                                    #    ###print('  Constraints added:', end=' ')
                                    #    addedFlag = True
                                    ####print('c'+str(c1.getName()), end=' ')
                        ###if addedFlag:
                            ###print()
                        
                        # add emerging constraints:
                        for i1 in [i for i in range(self.num_items) if i not in freq_const_id]:
                            c1 = self.constraints[i1+self.num_transactions+self.num_items]
                            if 1 in self.getDomain(i1).getValues():
                                #if 0 in self.getDomain(i1).getValues():
                                if c1 not in Q:
                                    Q.append(c1)
                                    #if not addedFlag:
                                    #    ###print('  Constraints added:', end=' ')
                                    #    addedFlag = True
                                    ####print('c'+str(c1.getName()), end=' ')
                        ###if addedFlag:
                            ###print()
                            
                        # add closure constraints:
                        for i1 in [i for i in range(self.num_items) if i not in freq_const_id]:
                            c1 = self.constraints[i1+self.num_transactions+2*self.num_items]
                            if 0 in self.getDomain(i1).getValues():
                                if c1 not in Q:
                                    Q.append(c1)
                                    #if not addedFlag:
                                    #    ###print('  Constraints added:', end=' ')
                                    #    addedFlag = True
                                    ####print('c'+str(c1.getName()), end=' ')
                        ###if addedFlag:
                            ###print()
                    else:
                        # add closure constraints:
                        for i1 in [i for i in range(self.num_items) if DB[c.getVariableRelativeIndex(),i] == 0]:
                            c1 = self.constraints[i1+self.num_transactions+2*self.num_items]
                            if 1 in self.getDomain(i1).getValues():
                                if c1 not in Q:
                                    Q.append(c1)
                                    #if not addedFlag:
                                    #    ###print('  Constraints added:', end=' ')
                                    #    addedFlag = True
                                    ####print('c'+str(c1.getName()), end=' ')
                        ###if addedFlag:
                            ###print()
                        
            
            # process frequency constraints:
            elif (c.getType() == 'frequency') or (c.getType() == 'emerging') or (c.getType() == 'closure'):
                if self.reviseUnary(self.getDomain(c.getVariableIndex()),c,0):
                    if self.getDomain(c.getVariableIndex()).isEmpty():
                        # print explanation:
                        #!self.failureExplanation(c.getFirstVariable())
                    
                        return False
                    
                    #addedFlag = False
                    if self.getDomain(c.getVariableIndex()).getValues()[0] != 1:                        
                        # add coverage constraints:
                        for i1 in [i for i in range(self.num_transactions) if DB2[i,c.getVariableRelativeIndex()] != 0]:
                            c1 = self.constraints[i1]
                            if 0 in self.getDomain(c1.getVariableIndex()).getValues():
                                if c1 not in Q:
                                    Q.append(c1)
                                    #if not addedFlag:
                                    #    ###print('  Constraints added:', end=' ')
                                    #    addedFlag = True
                                    ###print('c'+str(c1.getName()), end=' ')
                        ###if addedFlag:
                            ###print()
                            
                    else:
                        # add coverage constraints:
                        for i1 in [i for i in range(self.num_transactions) if DB2[i,c.getVariableRelativeIndex()] != 0]:
                            c1 = self.constraints[i1]
                            if 1 in self.getDomain(c1.getVariableIndex()).getValues():
                                if c1 not in Q:
                                    Q.append(c1)
                                    #if not addedFlag:
                                    #    ###print('  Constraints added:', end=' ')
                                    #    addedFlag = True
                                    ###print('c'+str(c1.getName()), end=' ')
                        ###if addedFlag:
                            ###print()
                            
                        # add size constraints:
                        for i1 in [i for i in range(self.num_transactions) if DB[i,c.getVariableRelativeIndex()] != 0]:
                            c1 = self.constraints[i1+self.num_transactions+3*self.num_items]
                            if 1 in self.getDomain(c1.getVariableIndex()).getValues():
                                if c1 not in Q:
                                    Q.append(c1)
                                    #if not addedFlag:
                                    #    ###print('  Constraints added:', end=' ')
                                    #    addedFlag = True
                                    ###print('c'+str(c1.getName()), end=' ')
                        ###if addedFlag:
                            ###print()

        return True
        
    def propagSuppress(self,Q,k):
        # process unary constraints:
        for c in Q:
            if c.getType() == 'unary':
                #!print('process unary constraint c'+str(c.getName()))
                if self.reviseUnary(self.getDomain(c.getVariableIndex()),c,k):
                    if self.getDomain(c.getVariableIndex()).isEmpty():
                        # print explanation:
                        #!self.failureExplanation(c.getFirstVariable())
                        
                        return False
                    
                    if c.getVariableType() == 'T':
                        if self.getDomain(c.getVariableIndex()).getValues()[0] != 1:
                            # add frequency constraints:
                            freq_const_id = [i for i in range(self.num_items) if DB[c.getVariableRelativeIndex(),i] != 0]
                            for i1 in freq_const_id:
                                c1 = self.constraints[i1+self.num_transactions]
                                if 1 in self.getDomain(i1).getValues():
                                    if c1 not in Q:
                                        Q.append(c1)
                                        #if not addedFlag:
                                        #    ###print('  Constraints added:', end=' ')
                                        #    addedFlag = True
                                        ####print('c'+str(c1.getName()), end=' ')
                            ###if addedFlag:
                                ###print()
                            
                            # add emerging constraints:
                            for i1 in [i for i in range(self.num_items) if i not in freq_const_id]:
                                c1 = self.constraints[i1+self.num_transactions+self.num_items]
                                if 1 in self.getDomain(i1).getValues():
                                    if c1 not in Q:
                                        Q.append(c1)
                                        #if not addedFlag:
                                        #    ###print('  Constraints added:', end=' ')
                                        #    addedFlag = True
                                        ####print('c'+str(c1.getName()), end=' ')
                            ###if addedFlag:
                                ###print()
                                
                            # add closure constraints:
                            for i1 in [i for i in range(self.num_items) if i not in freq_const_id]:
                                c1 = self.constraints[i1+self.num_transactions+2*self.num_items]
                                if 0 in self.getDomain(i1).getValues():
                                    if c1 not in Q:
                                        Q.append(c1)
                                        #if not addedFlag:
                                        #    ###print('  Constraints added:', end=' ')
                                        #    addedFlag = True
                                        ####print('c'+str(c1.getName()), end=' ')
                            ###if addedFlag:
                                ###print()
                        else:
                            # add closure constraints:
                            for i1 in [i for i in range(self.num_items) if DB[c.getVariableRelativeIndex(),i] == 0]:
                                c1 = self.constraints[i1+self.num_transactions+2*self.num_items]
                                if 1 in self.getDomain(i1).getValues():
                                    if c1 not in Q:
                                        Q.append(c1)
                                        #if not addedFlag:
                                        #    ###print('  Constraints added:', end=' ')
                                        #    addedFlag = True
                                        ####print('c'+str(c1.getName()), end=' ')
                            ###if addedFlag:
                                ###print()

                    elif c.getVariableType() == 'I':
                        #addedFlag = False
                        if self.getDomain(c.getVariableIndex()).getValues()[0] != 1:                      
                            # add coverage constraints:
                            for i1 in [i for i in range(self.num_transactions) if DB2[i,c.getVariableRelativeIndex()] != 0]:
                                c1 = self.constraints[i1]
                                if 0 in self.getDomain(c1.getVariableIndex()).getValues():
                                    if c1 not in Q:
                                        Q.append(c1)
                                        #if not addedFlag:
                                        #    ###print('  Constraints added:', end=' ')
                                        #    addedFlag = True
                                        ###print('c'+str(c1.getName()), end=' ')
                            ###if addedFlag:
                                ###print()
                        else:
                            # add coverage constraints:
                            for i1 in [i for i in range(self.num_transactions) if DB2[i,c.getVariableRelativeIndex()] != 0]:
                                c1 = self.constraints[i1]
                                if 1 in self.getDomain(c1.getVariableIndex()).getValues():
                                    if c1 not in Q:
                                        Q.append(c1)
                                        #if not addedFlag:
                                        #    ###print('  Constraints added:', end=' ')
                                        #    addedFlag = True
                                        ###print('c'+str(c1.getName()), end=' ')
                            ###if addedFlag:
                                ###print()
                                
                            # add size constraints:
                            for i1 in [i for i in range(self.num_transactions) if DB[i,c.getVariableRelativeIndex()] != 0]:
                                c1 = self.constraints[i1+self.num_transactions+3*self.num_items]
                                if 1 in self.getDomain(c1.getVariableIndex()).getValues():
                                    if c1 not in Q:
                                        Q.append(c1)
                                        #if not addedFlag:
                                        #    ###print('  Constraints added:', end=' ')
                                        #    addedFlag = True
                                        ###print('c'+str(c1.getName()), end=' ')
                            ###if addedFlag:
                                ###print()
                        
                Q.remove(c)
        
        while Q != []:
            c = Q.pop(0)
            
            # process coverage constraints:
            if (c.getType() == 'coverage') or (c.getType() == 'size'):
                #!print('process coverage constraint c'+str(c.getName()))
                if self.reviseUnary(self.getDomain(c.getVariableIndex()),c,k):
                    if self.getDomain(c.getVariableIndex()).isEmpty():
                        # print explanation:
                        #!self.failureExplanation(c.getFirstVariable())
                    
                        return False
                    
                    #addedFlag = False
                    if self.getDomain(c.getVariableIndex()).getValues()[0] != 1:
                        # add frequency constraints:
                        freq_const_id = [i for i in range(self.num_items) if DB[c.getVariableRelativeIndex(),i] != 0]
                        for i1 in freq_const_id:
                            c1 = self.constraints[i1+self.num_transactions]
                            if 1 in self.getDomain(i1).getValues():
                                if c1 not in Q:
                                    Q.append(c1)
                                    #if not addedFlag:
                                    #    ###print('  Constraints added:', end=' ')
                                    #    addedFlag = True
                                    ####print('c'+str(c1.getName()), end=' ')
                            ###if addedFlag:
                                ###print()
                        
                        # add emerging constraints:
                        for i1 in [i for i in range(self.num_items) if i not in freq_const_id]:
                            c1 = self.constraints[i1+self.num_transactions+self.num_items]
                            if 1 in self.getDomain(i1).getValues():
                                if c1 not in Q:
                                    Q.append(c1)
                                    #if not addedFlag:
                                    #    ###print('  Constraints added:', end=' ')
                                    #    addedFlag = True
                                    ####print('c'+str(c1.getName()), end=' ')
                        ###if addedFlag:
                            ###print()
                            
                        # add closure constraints:
                        for i1 in [i for i in range(self.num_items) if i not in freq_const_id]:
                            c1 = self.constraints[i1+self.num_transactions+2*self.num_items]
                            if 0 in self.getDomain(i1).getValues():
                                if c1 not in Q:
                                    Q.append(c1)
                                    #if not addedFlag:
                                    #    ###print('  Constraints added:', end=' ')
                                    #    addedFlag = True
                                    ####print('c'+str(c1.getName()), end=' ')
                        ###if addedFlag:
                            ###print()
                    else:
                        # add closure constraints:
                        for i1 in [i for i in range(self.num_items) if DB[c.getVariableRelativeIndex(),i] == 0]:
                            c1 = self.constraints[i1+self.num_transactions+2*self.num_items]
                            if 1 in self.getDomain(i1).getValues():
                                if c1 not in Q:
                                    Q.append(c1)
                                    #if not addedFlag:
                                    #    ###print('  Constraints added:', end=' ')
                                    #    addedFlag = True
                                    ####print('c'+str(c1.getName()), end=' ')
                        ###if addedFlag:
                            ###print()
                        
            
            # process frequency constraints:
            elif (c.getType() == 'frequency') or (c.getType() == 'emerging') or (c.getType() == 'closure'):
                #!print('process frequency constraint c'+str(c.getName()))
                if self.reviseUnary(self.getDomain(c.getVariableIndex()),c,k):
                    if self.getDomain(c.getVariableIndex()).isEmpty():
                        # print explanation:
                        #!self.failureExplanation(c.getFirstVariable())
                    
                        return False
                    
                    #addedFlag = False
                    if self.getDomain(c.getVariableIndex()).getValues()[0] != 1:                   
                        # add coverage constraints:
                        for i1 in [i for i in range(self.num_transactions) if DB2[i,c.getVariableRelativeIndex()] != 0]:
                            c1 = self.constraints[i1]
                            if 0 in self.getDomain(c1.getVariableIndex()).getValues():
                                if c1 not in Q:
                                    Q.append(c1)
                                    #if not addedFlag:
                                    #    ###print('  Constraints added:', end=' ')
                                    #    addedFlag = True
                                    ###print('c'+str(c1.getName()), end=' ')
                        ###if addedFlag:
                            ###print()
                    else:
                        # add coverage constraints:
                        for i1 in [i for i in range(self.num_transactions) if DB2[i,c.getVariableRelativeIndex()] != 0]:
                            c1 = self.constraints[i1]
                            if 1 in self.getDomain(c1.getVariableIndex()).getValues():
                                if c1 not in Q:
                                    Q.append(c1)
                                    #if not addedFlag:
                                    #    ###print('  Constraints added:', end=' ')
                                    #    addedFlag = True
                                    ###print('c'+str(c1.getName()), end=' ')
                        ###if addedFlag:
                            ###print()
                        
                        # add size constraints:
                        for i1 in [i for i in range(self.num_transactions) if DB[i,c.getVariableRelativeIndex()] != 0]:
                            c1 = self.constraints[i1+self.num_transactions+3*self.num_items]
                            if 1 in self.getDomain(c1.getVariableIndex()).getValues():
                                if c1 not in Q:
                                    Q.append(c1)
                                    #if not addedFlag:
                                    #    ###print('  Constraints added:', end=' ')
                                    #    addedFlag = True
                                    ###print('c'+str(c1.getName()), end=' ')
                        ###if addedFlag:
                            ###print()

        return True
    
    def restoreLevel(self,k):
        while self.removed_step[k] != []:
            v,b = self.removed_step[k].pop(0)
            self.restoreDomainValue(v,b)
        
    # search all solutions:
    def searchSolutionMAC(self,k):
        global numSolutions
        global stopFlag
        if (k==self.num_items) or stopFlag:
            return True
        else:
            consistent = True
            while consistent:
                v_i = self.not_instantiated.pop(0)
                a = self.getDomain(v_i.getIndex()).getValues()[0]

                newConstraintId = self.getMaxConstraintId()
                Q = [Constraint(newConstraintId,'unary',v_i,None,'=',a)]
                
                if self.propagSuppress(Q,k+1) and self.searchSolutionMAC(k+1):
                    if self.checkSolution():
                        numSolutions = numSolutions + 1

                self.restoreLevel(k+1)
                self.not_instantiated.append(v_i)

                newConstraintId = self.getMaxConstraintId()
                Q = [Constraint(newConstraintId,'unary',v_i,None,'!=',a)]

                if self.propagSuppress(Q,k):
                    consistent = True
                else:
                    consistent = False

        return consistent
        
    def failureExplanation(self,variableName):
        print('Expl(CSP->Failure) := {Domain '+variableName+'-> []}')
        self.printDomains()
            
    def successExplanation(self):
        print('Expl(CSP->Success) := {', end=' ')
        for i in range(len(self.variables)):
            if i < len(self.variables) - 1:
                print(self.variables[i].getName()+' = '+str(self.getDomain(i).getValues()), end=', ')
            else:
                print(self.variables[i].getName()+' = '+str(self.getDomain(i).getValues())+' }')
                
    def checkSuccess(self):
        check = True
        for i in range(len(self.domains)):
            if len(self.domains[i].getValues()) != 1:
                check = False
        if check:
            self.successExplanation()
        
    def supportExists(self,valueDx,domainY,constraint,orientation):
        exists = False
        for valueDy in domainY.getValues():
            if constraint.verify(valueDx,valueDy,orientation):
                exists = True
            
        return exists
       
    def reviseUnary(self,domainX,constraint,k):
        eliminated = False
        domainValues = domainX.getValues()[:]
        for valueDx in domainValues:
            if not constraint.verify(valueDx,None,'forward',self.items_min,self.items_max,self.transactions_min,self.transactions_max):
                self.removeDomainValue(constraint.getFirstVariable(),valueDx)
                if k > 0:
                    self.removed_step[k].append((constraint.getFirstVariable(),valueDx))
                eliminated = True
                
        return eliminated

    def reviseOriented(self,domainX,domainY,constraint,orientation,k):     
        # initialization:
        eliminated = False
        
        # Explanations for the values which are still present:
        domainValues = domainX.getValues()[:]
        for valueDx in domainValues:
            if not self.supportExists(valueDx,domainY,constraint,orientation):
                domainX.removeValue(valueDx)
                if k > 0:
                    self.removed_step[k].append((self.getVariableName(domainX),valueDx))
                    
                eliminated = True
                
        return eliminated
                           
    def printAllExplanations(self):
        for variableName,value in self.all_explanations:
            allExplanations = self.all_explanations[variableName,value][:]
            
            print('AllExpl('+variableName+'≠'+str(value)+') := {', end=' ')
            
            n = len(allExplanations)
            for i in range(n):
                m = len(allExplanations[i])
                for j in range(m):
                    if (j == 0) and (m != 1):
                        print('{c'+str(allExplanations[i][j])+',', end = ' ')
                    elif (j == 0) and (m == 1):
                        print('{c'+str(allExplanations[i][j])+'}', end = ' ')
                    elif (j != 0) and (j != m-1):
                        print('c'+str(allExplanations[i][j])+',', end = ' ')
                    elif (j != 0) and (j == m-1):
                        print('c'+str(allExplanations[i][j])+'}', end = ' ')
                if i != n-1:
                    print(',', end = ' ')
                else:
                    print('}')

def loadDBFromFile(file):
    f = open('data/'+file, 'r')
   
    numTrans = 0
    maxItem = 0
    for line in f:
        if line !='\n' and not '@' in line:
            data = line.strip().split(' ')
            maxElement = np.max([int(i) for i in data[:-1]])
            if maxElement > maxItem:
                maxItem = maxElement
            numTrans += 1
    f.close()
    
    numItems = maxItem + 1
    
    dataBase = np.zeros((numTrans,numItems),int)
    dataBase0 = np.zeros((numTrans,numItems),int)
    dataBase1 = np.zeros((numTrans,numItems),int)
    
    f = open('data/'+file, 'r')
   
    count = 0
    count_0 = 0
    count_1 = 0
    for line in f:
        if line !='\n' and not '@' in line:
            data = line.strip().split(' ')
            if data[-1] == '1':
                for item in data[:-1]:
                    dataBase1[count][int(item)] = True
                count_1+=1
            elif data[-1] == '0':
                for item in data[:-1]:
                    dataBase0[count][int(item)] = True
                count_0+=1
            for item in data[:-1]:
                    dataBase[count][int(item)] = True
            count += 1
    f.close()
    
    return dataBase,dataBase0,dataBase1,count_0,count_1

def loadDBsFromFiles(fileNameTr,fileNameAct):
    # parsing molecule activity file:
    f = open('data/'+fileNameAct, 'r')
    
    line = f.readline()
    numTrans = 0
    moleculeActivity = []
    for line in f:
        if line !='\n':
            data = line.strip().split(';')
            moleculeActivity.append(int(data[1]))
            numTrans += 1
    f.close()
    
    # parsing pharmacophore file:
    f = open('data/'+fileNameTr, 'r')
    
    line = f.readline()
    data = line.strip('\n').strip(';').split(';')
    
    pharma_dict = {}
    count = 0
    for item in data:
        pharma_dict[count] = ''.join(item.split(' '))
        count+=1
    
    numItems = len(data)
    dataBase = np.zeros((numTrans,numItems),int)
    dataBase0 = np.zeros((numTrans,numItems),int)
    dataBase1 = np.zeros((numTrans,numItems),int)
    
    count_0 = 0
    count_1 = 0
    for line in f:
        if line !='\n':
            data = line.strip('\n').strip(';').split(';')
            if moleculeActivity[int(data[0])] == 0:
                dataBase0[int(data[0])][:] = [item=='1' for item in data[1:]]
                count_0+=1
            else:
                dataBase1[int(data[0])][:] = [item=='1' for item in data[1:]]
                count_1+=1
            dataBase[int(data[0])][:] = [item=='1' for item in data[1:]]
    f.close()
    
    return dataBase,dataBase0,dataBase1,count_0,count_1,np.array(moleculeActivity),pharma_dict

# load DB:
DB,DB0,DB1,DN,DP,DC,PHARMA = loadDBsFromFiles(fileNameTr,fileNameAct)
#DB,DB0,DB1,DN,DP, = loadDBFromFile(fileName)
DB2 = 1-DB
CD = 1-DC

# determine the numbers of items and transactions:
numItems = np.shape(DB)[1]
numTrans = np.shape(DB)[0]

DPD = DP/float(numTrans)
DND = DN/float(numTrans)
D = numTrans

print('numItems: ',numItems)
print('numTrans: ',numTrans)

print('THETA: ',THETA)
print('THETA2: ',THETA2)
print('THETA3: ',THETA3)

# load variables and domains:
variables = []
domains = []
constraints = []

# add variables:
# items:
for i in range(numItems):
    variables.append(Variable('I'+str(i+1),i,i))

# transactions:
for i in range(numTrans):
    variables.append(Variable('T'+str(i+1),i,i+numItems))

# add domains:
for i in range(numItems+numTrans):
    domains.append(Domain([1,0]))

# add coverage constraints:
for i in range(numTrans):
    constraints.append(Constraint(i+1,'coverage',variables[i+numItems]))

# add frequency constraints:
for i in range(numItems):
    constraints.append(Constraint(numTrans+i+1,'frequency',variables[i]))

# add emerging constraints:
for i in range(numItems):
    constraints.append(Constraint(numTrans+numItems+i+1,'emerging',variables[i]))

# add closure constraints:
for i in range(numItems):
    constraints.append(Constraint(numTrans+2*numItems+i+1,'closure',variables[i]))
    
# add size constraints:
for i in range(numTrans):
    constraints.append(Constraint(numTrans+3*numItems+i+1,'size',variables[i+numItems]))

# load the CSP:
task = CSP(numItems,numTrans,variables,domains,constraints)

f_out = open('results/'+fileName+'_emerging_closed_'+str(THETA)+'_p'+str(THETA2)+'_f'+str(THETA3)+'.csv','w')

f_out.write('N,Motif,Molécules actives,Molécules inactives,Chi-square\n')

start = timeit.default_timer() # Initialize timer to compute the running time

numSolutions = 0
stopFlag = False

constr = []

# construct the search tree:
result = task.mac()

if result:
    task.checkSuccess()
    
stop = timeit.default_timer() # stop time counting
print('Number of solutions: ',numSolutions)
print('Positive solutions: ',task.num_pos_patterns)
print('Negative solutions: ',task.num_neg_patterns)
print('Number of constraints: ',task.getMaxConstraintId()-1)
print('Iteration time, sec: '+str(np.round(stop-start,2)))
f_out.close()
