class GenomicRegion:

  def __init__(self,chr,initial,final,name=None,orientation=None,data=None):
    self.chr=chr
    self.initial=initial
    self.final=final
    self.name=name
    self.orientation=orientation
    self.data=data

  def __len__(self):
    return self.final-self.initial


  def __str__(self):
    s= self.chr+'\t'+str(self.initial)+'\t'+str(self.final)
    if self.name != None:
      s=s+'\t'+self.name
    if self.orientation != None:
      s=s+'\t'+self.orientation
    if self.data != None:
      s=s+'\t'+self.data
    return s

  def extend(self,left,right):
    self.initial-=left
    self.final+=right    

  def overlapp(self,region):
    if region.chr == self.chr:
        if self.initial < region.initial:
            if self.final >= region.initial:
                return True
        else:
            if self.initial <= region.final:
                return True
    return False

                    
  def __repr__(self):
    return self.chr+','+str(self.initial)+','+str(self.final)
  

  def __cmp__(x,y):
    '''  The return value is negative if x < y, zero if x == y and strictly positive if x > y. '''
    if x.chr<y.chr:
      return -1
    elif x.chr>y.chr:
      return 1
    else:
      if x.initial < y.initial:
        return -1
      elif x.initial > y.final:
        return 1
      else:
        if x.final < y.final:
          return -1
        elif x.final > y.final:
          return 1
        else:
          return 0



