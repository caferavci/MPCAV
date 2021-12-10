# -*- coding: utf-8 -*-
"""
Created on Tue Aug  3 02:52:31 2021

@authors: 
Claudio Roncoli and Cafer Avci
"""

import sys
import mosek
import numpy as np

# Since the value of infinity is ignored, we define it solely
# for symbolic purposes
inf = 0.0

# Define a stream printer to grab output from MOSEK
def streamprinter(text):
	sys.stdout.write(text)
	sys.stdout.flush()

class Cell:
	def __init__(self, id, queue):
		self.id = id
		self.enteringLinks = []
		self.exitingLinks = []
		self.queue = queue
		self.h = 0
		self.divergingIndex = 0
		
	def addEnteringLink(self, link):
		link.downstreamCell = self
		self.enteringLinks.append(link)
	
	def addExitingLink(self, link):
		link.upstreamCell = self
		self.exitingLinks.append(link)
	
	def setStatus(self):
		self.status = 0
		
		if len(self.exitingLinks) > 1:
			self.status = 1
		if len(self.enteringLinks) > 1:
			self.status = 2
		if len(self.enteringLinks) == 0:
			self.status = 3
		if len(self.exitingLinks) == 0:
			self.status = 4
			
	def setPenaltyLabel(self, h):
		self.h = h
		
	def setPenaltyRecursive(self):
		for el in self.enteringLinks:
			el.upstreamCell.setPenaltyLabel(self.h+1)
			el.upstreamCell.setPenaltyRecursive()
			
	def setDivergingIndex(self, divergingIndex):
		self.divergingIndex = divergingIndex

class Link:
	def __init__(self, id, type):
		self.id = id
		self.type = type # 0: ordinary, 1: exit, 2: entering demand
	
	def setUpstreamCell(self, cell):
		self.upstreamCell = cell

	def setDownstreamCell(self, cell):
		self.downstreamCell = cell
		
def mainOptimisation(id):

	numCells = id.numCells
	numLinks = id.numLinks
	numSources = id.numSources
	
	linkPattern = id.linkPattern
	
	vMax = id.vMax
	minHeadway = id.minHeadway
	maxQueue = id.maxQueue
	
	horizon = id.horizon
	
	T = id.T
	D = id.D
	L = id.L
	
	weightQ = id.weightQ
	
	rhoInit = id.rhoInit
	queueInit = id.queueInit
	
	demand = id.demand
		
	maxH = 3600/np.array(minHeadway)
	
	cdf = id.capacityDropFactor	
	
	cells = []
	
	# create all the cells
	for i in range(numCells):
		cell = Cell(i,0)
		cells.append(cell)
	
	# create links and add them to cells
	for i in range(len(linkPattern)):
		lp = linkPattern[i]

		l = Link(i,0)
		
		cells[lp[0]].addExitingLink(l)
		cells[lp[1]].addEnteringLink(l)
	
	# define cell status: sink = 4, source = 3, merging = 2, diverging = 1, ordinary = 0
	cellStatus = [0] * 5
	
	numStates = numCells
	numInputs = numLinks
	
	indDiverging = 0
	
	for i in range(len(cells)):
		c = cells[i]
		
		c.setStatus()
		
		# add exiting link at sink cells
		if c.status == 4:
			l = Link(numInputs,1)
			c.addExitingLink(l)
			
			numInputs += 1
			
			sinkId = c.id
			
		# add entering link and queue at source cells
		if c.status == 3:
			l = Link(numInputs,2)
			c.addEnteringLink(l)
			
			numInputs += 1
			
			cell = Cell(numStates,1)
			cell.addExitingLink(l)
			cells.append(cell)
			numStates += 1
		
		# add exiting link at sink cells
		if c.status == 1:
			c.setDivergingIndex(indDiverging)
			indDiverging += 1
		
		cellStatus[c.status] += 1
	
	print('Source cells:', str(cellStatus[3]), '\nSink cells:', str(cellStatus[4]), '\nDiverging cells:', str(cellStatus[1]), '\nMerging cells:', str(cellStatus[2]), '\nOrdinary cells:', str(cellStatus[0]))
	
	blockSize = numInputs + numStates
	
	numvar = blockSize * horizon
	
	if id.avoidHB:
		numvar += 1
		
		cells[sinkId].setPenaltyRecursive()			
		hbsub = []
		hbval = []
	
	numcon = 0
	
	# constraints matrices
	asub = []
	aval = []
	
	bkc = []
	blc = []
	buc = []
	
	# bounds arrays
	bkx = [0] * numvar
	blx = [0] * numvar
	bux = [0] * numvar
	
	# cost matrices
	linCost = [0] * numvar
	qsubi = []
	qsubj = []
	qval = []
	
	for k in range(horizon):
		offset = k * blockSize
		
		for i in range(numCells):
			cell = cells[i]
			
			if k == 0:
				if cell.status == 4 or cell.status == 3 or cell.status == 0:
					# density dynamics
					asub.append([numInputs+cell.id, cell.enteringLinks[0].id, cell.exitingLinks[0].id])
					aval.append([1, -T/D, T/D])

					bkc.append(mosek.boundkey.fx)
					blc.append(rhoInit[cell.id])
					buc.append(rhoInit[cell.id])
					
					numcon += 1
				
				if cell.status == 1:
					# density dynamics
					asub.append([numInputs+cell.id, cell.enteringLinks[0].id, cell.exitingLinks[0].id, cell.exitingLinks[1].id])
					aval.append([1, -T/D, T/D, T/D])
					
					bkc.append(mosek.boundkey.fx)
					blc.append(0)
					buc.append(0)
					
					numcon += 1
					
					# demand undercritical
					asub.append([offset+cell.exitingLinks[0].id, offset+cell.exitingLinks[1].id])
					aval.append([1, 1])
					
					bkc.append(mosek.boundkey.up)
					blc.append(-inf)
					buc.append(vMax[cell.id] * rhoInit[cell.id])
					
					numcon += 1
					
					# demand overcritical
					asub.append([offset+cell.exitingLinks[0].id, offset+cell.exitingLinks[1].id])
					aval.append([1, 1])
					
					bkc.append(mosek.boundkey.up)
					blc.append(-inf)
					buc.append( ( vMax[cell.id] * maxH[cell.id] ) / ( vMax[cell.id] + L * maxH[cell.id] ) )
					
					numcon += 1
					
					for l in cell.exitingLinks:
						# supply overcritical
						asub.append([offset+l.id, offset+numInputs+l.downstreamCell.id])
						aval.append([1, L*maxH[l.downstreamCell.id]])
					
						bkc.append(mosek.boundkey.up)
						blc.append(-inf)
						buc.append( maxH[cell.id] )
					
						numcon += 1
						
						# supply undercritical and overcritical (as bound)
						bkx[offset+l.id] = mosek.boundkey.ra
						blx[offset+l.id] = 0.0
						bux[offset+l.id] = min ( ( vMax[l.downstreamCell.id] * maxH[l.downstreamCell.id] ) / ( vMax[l.downstreamCell.id] + L * maxH[l.downstreamCell.id] ), -L*maxH[l.downstreamCell.id]*rhoInit[l.downstreamCell.id]+maxH[l.downstreamCell.id])					
					
					# turning rates: works ONLY if there are two outbound links!
					if id.fixedTurning:
						asub.append([offset+cell.exitingLinks[0].id, offset+cell.exitingLinks[1].id])
						aval.append([1-id.turning[cell.divergingIndex,k], -id.turning[cell.divergingIndex,k]])
					
						bkc.append(mosek.boundkey.fx)
						blc.append(0)
						buc.append(0)
					
						numcon += 1
					
				if cell.status == 2:
					# density dynamics
					asub.append([numInputs+cell.id, cell.enteringLinks[0].id, cell.enteringLinks[1].id, cell.exitingLinks[0].id])
					aval.append([1, -T/D, -T/D, T/D])
					
					bkc.append(mosek.boundkey.fx)
					blc.append(0)
					buc.append(0)
					
					numcon += 1
					
					# supply undercritical and overcritical for entering flows & capacity drop (outflow of the current cell as function of densities in entering cells)
					upCell = cell.enteringLinks[0].upstreamCell
					capDrop0 = cdf*vMax[cell.id]*maxH[cell.id] / (vMax[cell.id]+L*maxH[cell.id]) + (1-cdf)*vMax[cell.id]*maxH[cell.id]*(vMax[upCell.id]+L*maxH[upCell.id]) / ( vMax[upCell.id]*(vMax[cell.id]+L*maxH[cell.id]) ) - ( (1-cdf)*vMax[cell.id]*L*maxH[cell.id]*(vMax[upCell.id]+L*maxH[upCell.id]) ) / ( vMax[upCell.id] * (vMax[cell.id]+L*maxH[cell.id]) ) * rhoInit[upCell.id]
					
					upCell = cell.enteringLinks[1].upstreamCell
					capDrop1 = cdf*vMax[cell.id]*maxH[cell.id] / (vMax[cell.id]+L*maxH[cell.id]) + (1-cdf)*vMax[cell.id]*maxH[cell.id]*(vMax[upCell.id]+L*maxH[upCell.id]) / ( vMax[upCell.id]*(vMax[cell.id]+L*maxH[cell.id]) ) - ( (1-cdf)*vMax[cell.id]*L*maxH[cell.id]*(vMax[upCell.id]+L*maxH[upCell.id]) ) / ( vMax[upCell.id] * (vMax[cell.id]+L*maxH[cell.id]) ) * rhoInit[upCell.id]
					
					asub.append([cell.enteringLinks[0].id, cell.enteringLinks[1].id])
					aval.append([1, 1])
					
					bkc.append(mosek.boundkey.up)
					blc.append(-inf)
					buc.append( min( -L*maxH[cell.id]*rhoInit[cell.id]+maxH[cell.id], ( vMax[cell.id] * maxH[cell.id] ) / ( vMax[cell.id] + L * maxH[cell.id] ), capDrop0, capDrop1 ) )
					
					numcon += 1			
					
					
				if cell. status == 3 or cell.status == 0:
					downstreamCell = cell.exitingLinks[0].downstreamCell
					# demand and supply, undercritical and overcritical (as bound)
					bkx[cell.exitingLinks[0].id] = mosek.boundkey.ra
					blx[cell.exitingLinks[0].id] = 0.0
					bux[cell.exitingLinks[0].id] = min( vMax[cell.id] * rhoInit[cell.id], ( vMax[cell.id] * maxH[cell.id] ) / ( vMax[cell.id] + L * maxH[cell.id] ), -L*maxH[downstreamCell.id]*rhoInit[downstreamCell.id]+maxH[downstreamCell.id], ( vMax[downstreamCell.id] * maxH[downstreamCell.id] ) / ( vMax[downstreamCell.id] + L * maxH[downstreamCell.id] ) )
				
				if cell.status == 3:
					queueCell = cell.enteringLinks[0].upstreamCell
					# queue dynamics
					asub.append([numInputs+queueCell.id, cell.enteringLinks[0].id])
					aval.append([1, T])
					
					bkc.append(mosek.boundkey.fx)
					blc.append(queueInit[queueCell.id-numCells] + T*demand[queueCell.id-numCells,k])
					buc.append(queueInit[queueCell.id-numCells] + T*demand[queueCell.id-numCells,k])
					
					numcon += 1
					
					# supply undercritical and overcritical for entering flow (as bound)
					bkx[cell.enteringLinks[0].id] = mosek.boundkey.ra
					blx[cell.enteringLinks[0].id] = 0.0
					bux[cell.enteringLinks[0].id] = min ( -L*maxH[cell.id]*rhoInit[cell.id]+maxH[cell.id], ( vMax[cell.id] * maxH[cell.id] ) / ( vMax[cell.id] + L * maxH[cell.id] ) )
				
				if cell.status == 4 or cell.status == 2:
					# demand undercritical and overcritical (as bound)
					bkx[cell.exitingLinks[0].id] = mosek.boundkey.ra
					blx[cell.exitingLinks[0].id] = 0.0
					bux[cell.exitingLinks[0].id] = min( vMax[cell.id] * rhoInit[cell.id], ( vMax[cell.id] * maxH[cell.id] ) / ( vMax[cell.id] + L * maxH[cell.id] ) )
				
				# quadratic cost
				for l in cell.exitingLinks:
					qsubi.append(l.id)
					qsubj.append(l.id)
					qval.append(weightQ*2)	

				if cell.status == 3:
					qsubi.append(cell.enteringLinks[0].id)
					qsubj.append(cell.enteringLinks[0].id)
					qval.append(weightQ*2)
					
			else:
				if cell.status == 4 or cell.status == 3 or cell.status == 0:
					# density dynamics
					asub.append([offset+numInputs+cell.id, offset+numInputs+cell.id-blockSize, offset+cell.enteringLinks[0].id, offset+cell.exitingLinks[0].id])
					aval.append([1, -1, -T/D, T/D])
					
					bkc.append(mosek.boundkey.fx)
					blc.append(0)
					buc.append(0)
					
					numcon += 1
				
				if cell.status == 1:
					# density dynamics
					asub.append([offset+numInputs+cell.id, offset+numInputs+cell.id-blockSize, offset+cell.enteringLinks[0].id, offset+cell.exitingLinks[0].id, offset+cell.exitingLinks[1].id])
					aval.append([1, -1, -T/D, T/D, T/D])
					
					bkc.append(mosek.boundkey.fx)
					blc.append(0)
					buc.append(0)
					
					numcon += 1
					
					# demand undercritical
					asub.append([offset+cell.exitingLinks[0].id, offset+cell.exitingLinks[1].id, offset+numInputs+cell.id-blockSize])
					aval.append([1, 1, -vMax[cell.id]])
					
					bkc.append(mosek.boundkey.up)
					blc.append(-inf)
					buc.append(0)
					
					numcon += 1
					
					# demand overcritical
					asub.append([offset+cell.exitingLinks[0].id, offset+cell.exitingLinks[1].id])
					aval.append([1, 1])
					
					bkc.append(mosek.boundkey.up)
					blc.append(-inf)
					buc.append( ( vMax[cell.id] * maxH[cell.id] ) / ( vMax[cell.id] + L * maxH[cell.id] ) )
					
					numcon += 1
					
					for l in cell.exitingLinks:
						# supply overcritical
						asub.append([offset+l.id, offset+numInputs+l.downstreamCell.id])
						aval.append([1, L*maxH[l.downstreamCell.id]])
					
						bkc.append(mosek.boundkey.up)
						blc.append(-inf)
						buc.append( maxH[cell.id] )
					
						numcon += 1
						
						# supply undercritical (as bound)
						bkx[offset+l.id] = mosek.boundkey.ra
						blx[offset+l.id] = 0.0
						bux[offset+l.id] = ( vMax[l.downstreamCell.id] * maxH[l.downstreamCell.id] ) / ( vMax[l.downstreamCell.id] + L * maxH[l.downstreamCell.id] )					
					
					# turning rates: works ONLY if there are two outbound links!
					if id.fixedTurning:
						asub.append([offset+cell.exitingLinks[0].id, offset+cell.exitingLinks[1].id])
						aval.append([1-id.turning[cell.divergingIndex,k], -id.turning[cell.divergingIndex,k]])
					
						bkc.append(mosek.boundkey.fx)
						blc.append(0)
						buc.append(0)
					
						numcon += 1
						
				if cell.status == 2:
					# density dynamics
					asub.append([offset+numInputs+cell.id, offset+numInputs+cell.id-blockSize, offset+cell.enteringLinks[0].id, offset+cell.enteringLinks[1].id, offset+cell.exitingLinks[0].id])
					aval.append([1, -1, -T/D, -T/D, T/D])
					
					bkc.append(mosek.boundkey.fx)
					blc.append(0)
					buc.append(0)
					
					numcon += 1
					
					# supply overcritical for entering flows
					asub.append([offset+cell.enteringLinks[0].id, offset+cell.enteringLinks[1].id, offset+numInputs+cell.id])
					aval.append([1, 1, L*maxH[cell.id]])
					
					bkc.append(mosek.boundkey.up)
					blc.append(-inf)
					buc.append(maxH[cell.id])
					
					numcon += 1
					
					# supply undercritical for entering flows
					asub.append([offset+cell.enteringLinks[0].id, offset+cell.enteringLinks[1].id])
					aval.append([1, 1])
					
					bkc.append(mosek.boundkey.up)
					blc.append(-inf)
					buc.append( ( vMax[cell.id] * maxH[cell.id] ) / ( vMax[cell.id] + L * maxH[cell.id] ) )
					
					numcon += 1
					
					# capacity drop (outflow of the current cell as function of densities in entering cells)
					upCell = cell.enteringLinks[0].upstreamCell
					asub.append([offset+cell.enteringLinks[0].id, offset+cell.enteringLinks[1].id, offset+numInputs+upCell.id])
					aval.append([1, 1, ( (1-cdf)*vMax[cell.id]*L*maxH[cell.id]*(vMax[upCell.id]+L*maxH[upCell.id]) ) / ( vMax[upCell.id] * (vMax[cell.id]+L*maxH[cell.id]) )])
					
					bkc.append(mosek.boundkey.up)
					blc.append(-inf)
					buc.append( cdf*vMax[cell.id]*maxH[cell.id] / (vMax[cell.id]+L*maxH[cell.id]) + (1-cdf)*vMax[cell.id]*maxH[cell.id]*(vMax[upCell.id]+L*maxH[upCell.id]) / ( vMax[upCell.id]*(vMax[cell.id]+L*maxH[cell.id]) ) )
					
					numcon += 1
					
					upCell = cell.enteringLinks[1].upstreamCell
					asub.append([offset+cell.enteringLinks[0].id, offset+cell.enteringLinks[1].id, offset+numInputs+upCell.id])
					aval.append([1, 1, ( (1-cdf)*vMax[cell.id]*L*maxH[cell.id]*(vMax[upCell.id]+L*maxH[upCell.id]) ) / ( vMax[upCell.id] * (vMax[cell.id]+L*maxH[cell.id]) )])
					
					bkc.append(mosek.boundkey.up)
					blc.append(-inf)
					buc.append( cdf*vMax[cell.id]*maxH[cell.id] / (vMax[cell.id]+L*maxH[cell.id]) + (1-cdf)*vMax[cell.id]*maxH[cell.id]*(vMax[upCell.id]+L*maxH[upCell.id]) / ( vMax[upCell.id]*(vMax[cell.id]+L*maxH[cell.id]) ) )
					
					numcon += 1
					
					
				if cell.status == 4 or cell.status == 3 or cell.status == 2 or cell.status == 0:
					# demand undercritical
					asub.append([offset+cell.exitingLinks[0].id, offset+numInputs+cell.id-blockSize])
					aval.append([1, -vMax[cell.id]])
					
					bkc.append(mosek.boundkey.up)
					blc.append(-inf)
					buc.append(0)
					
					numcon += 1
												
				if cell.status == 3 or cell.status == 0:
					downstreamCell = cell.exitingLinks[0].downstreamCell
					#supply overcritical
					asub.append([offset+cell.exitingLinks[0].id, offset+numInputs+downstreamCell.id])
					aval.append([1, L*maxH[downstreamCell.id]])
					
					bkc.append(mosek.boundkey.up)
					blc.append(-inf)
					buc.append(maxH[downstreamCell.id])
					
					numcon += 1
				
					# demand overcritical and supply undercritical (as bound)
					bkx[offset+cell.exitingLinks[0].id] = mosek.boundkey.ra
					blx[offset+cell.exitingLinks[0].id] = 0.0
					bux[offset+cell.exitingLinks[0].id] = min( ( vMax[cell.id] * maxH[cell.id] ) / ( vMax[cell.id] + L * maxH[cell.id] ), ( vMax[downstreamCell.id] * maxH[downstreamCell.id] ) / ( vMax[downstreamCell.id] + L * maxH[downstreamCell.id] ) )
					
					# bound on headway active on the demand term. ACTIVATE IF HEADWAY EXCEEDS THE LOWER BOUND
					asub.append([offset+cell.exitingLinks[0].id, offset+numInputs+cell.id])
					aval.append([1, L*maxH[cell.id]])
					
					bkc.append(mosek.boundkey.up)
					blc.append(-inf)
					buc.append(maxH[cell.id])
					
					numcon += 1
					
				if cell.status == 3:
					queueCell = cell.enteringLinks[0].upstreamCell
					# queue dynamics
					asub.append([offset+numInputs+queueCell.id, offset+numInputs+queueCell.id-blockSize, offset+cell.enteringLinks[0].id])
					aval.append([1, -1, T])
					
					bkc.append(mosek.boundkey.fx)
					blc.append(T*demand[queueCell.id-numCells,k])
					buc.append(T*demand[queueCell.id-numCells,k])
					
					numcon += 1
					
					# supply overcritical for entering flow
					asub.append([offset+cell.enteringLinks[0].id, offset+numInputs+cell.id])
					aval.append([1, L*maxH[cell.id]])
					
					bkc.append(mosek.boundkey.up)
					blc.append(-inf)
					buc.append(maxH[cell.id])
					
					numcon += 1
					
					# supply undercritical for entering flow (as bound)
					bkx[offset+cell.enteringLinks[0].id] = mosek.boundkey.ra
					blx[offset+cell.enteringLinks[0].id] = 0.0
					bux[offset+cell.enteringLinks[0].id] = ( vMax[cell.id] * maxH[cell.id] ) / ( vMax[cell.id] + L * maxH[cell.id] )
				
				if cell.status == 4 or cell.status == 2:
					# demand overcritical (as bound)
					bkx[offset+cell.exitingLinks[0].id] = mosek.boundkey.ra
					blx[offset+cell.exitingLinks[0].id] = 0.0
					bux[offset+cell.exitingLinks[0].id] = ( vMax[cell.id] * maxH[cell.id] ) / ( vMax[cell.id] + L * maxH[cell.id] )
					
				# quadratic cost
				for l in cell.exitingLinks:
					qsubi.append(offset+l.id)
					qsubj.append(offset+l.id)
					if k < horizon-1:
						qval.append(weightQ*4)
					else:
						qval.append(weightQ*2)
					
					qsubi.append(offset+l.id)
					qsubj.append(offset+l.id-blockSize)
					qval.append(weightQ*-2)

				if cell.status == 3:
					qsubi.append(offset+cell.enteringLinks[0].id)
					qsubj.append(offset+cell.enteringLinks[0].id)
					if k < horizon-1:
						qval.append(weightQ*4)
					else:
						qval.append(weightQ*2)
					
					qsubi.append(offset+cell.enteringLinks[0].id)
					qsubj.append(offset+cell.enteringLinks[0].id-blockSize)
					qval.append(-weightQ*2)					
					
					
			# all steps
			#density bounds
			bkx[offset+numInputs+cell.id] = mosek.boundkey.lo
			blx[offset+numInputs+cell.id] = 0.0
			bux[offset+numInputs+cell.id] = inf
			
			linCost[offset+numInputs+cell.id] = T*D
			
			# queue bounds			
			if cell.status == 3:
				queueCell = cell.enteringLinks[0].upstreamCell
				
				bkx[offset+numInputs+queueCell.id] = mosek.boundkey.ra
				blx[offset+numInputs+queueCell.id] = 0.0
				bux[offset+numInputs+queueCell.id] = maxQueue[queueCell.id-numCells]
				
				linCost[offset+numInputs+queueCell.id] = T
			
			if id.avoidHB and cell.status != 4:
				hbsub.append(offset+numInputs+cell.id)
				hbval.append(cell.h)
			
	if id.avoidHB:
		# linear constraint
		hbsub.append(numvar-1)
		hbval.append(-id.largeM)

		asub.append(hbsub)
		aval.append(hbval)		
		
		bkc.append(mosek.boundkey.up)
		blc.append(-inf)
		buc.append(0)
					
		numcon += 1
					
		# constraint and cost for m
		bkx[numvar-1] = mosek.boundkey.up
		blx[numvar-1] = -inf
		bux[numvar-1] = id.smallM
				
		linCost[numvar-1] = 1
	
	# Make mosek environment
	with mosek.Env() as env:
		# Create a task object
		with env.Task(0, 0) as task:
			# Attach a log stream printer to the task
			task.set_Stream(mosek.streamtype.log, streamprinter)
			
			# Append 'numcon' empty constraints.
			# The constraints will initially have no bounds.
			task.appendcons(numcon)

			# Append 'numvar' variables.
			# The variables will initially be fixed at zero (x=0).
			task.appendvars(numvar)

			for j in range(numvar):
				task.putvarbound(j, bkx[j], blx[j], bux[j])
				# Set the linear term c_j in the objective.
				task.putcj(j, 1000*linCost[j])
				
			for j in range(numcon):
				task.putarow(j,				  # Variable (row) index.
							 asub[j],			# column index of non-zeros in column j.
							 aval[j])			# Non-zero Values of column j.
				# Set the bounds on constraints
				task.putconbound(j, bkc[j], blc[j], buc[j])
			
			# quadratic cost
			task.putqobj(qsubi, qsubj, qval)
			
			# Input the objective sense (minimize/maximize)
			task.putobjsense(mosek.objsense.minimize)

			# Solve the problem
			task.optimize()
			# Print a summary containing information
			# about the solution for debugging purposes
			task.solutionsummary(mosek.streamtype.msg)

			# Get status information about the solution
			solsta = task.getsolsta(mosek.soltype.itr)

			if (solsta == mosek.solsta.optimal or solsta == mosek.solsta.near_optimal):
				xx = [0.] * numvar
				task.getxx(mosek.soltype.itr, xx)
				
				# for i in range(numvar):
					# print("x[" + str(i) + "]=" + str(xx[i]))
					
				rho = np.zeros((numCells,horizon+1))
				qu = np.zeros((cellStatus[3],horizon+1))
				q = np.zeros((numLinks+cellStatus[4], horizon))
				r = np.zeros((cellStatus[3],horizon))
				speed = np.zeros((numCells,horizon))
				headway = np.zeros((numCells,horizon))
				ratio = np.ones((numLinks+cellStatus[4], horizon))
				ttInst = np.zeros((len(id.paths),horizon))
				
				tts = 0
				
				j = 0
				for i in range(numCells):
					cell = cells[i]
					
					rho[i,0] = rhoInit[i]
					
					tts += T*D*rho[i,0]
					
					if cell.status == 3:
						qu[j,0] = queueInit[j]
						
						tts += T*qu[j,0]
						
						j += 1
				
				for k in range(horizon):
					offset = k * blockSize
					
					j = 0
					
					for i in range(numCells):
						cell = cells[i]
						
						rho[i,k+1] = xx[offset+numInputs+i]
					
						tts += T*D*rho[i,k+1]
												
						if cell.status == 3:
							queueCell = cell.enteringLinks[0].upstreamCell
							
							qu[j,k+1] = xx[offset+numInputs+queueCell.id]
							
							tts += T*qu[j,k+1]
							
							r[j,k] = xx[offset+cell.enteringLinks[0].id]
							
							j += 1
						
						totFlow = 0
						if cell.status != 4:
							for l in cell.exitingLinks:
								q[l.id,k] = xx[offset+l.id]
								totFlow += q[l.id,k]
						else:
							l = cell.exitingLinks[0]
							q[numLinks,k] = xx[offset+l.id]
							totFlow += q[numLinks,k]
						
						if cell.status == 1:
							for l in cell.exitingLinks:
								ratio[l.id,k] = q[l.id,k] / totFlow
						
						if rho[i,k] > 0.1:
							speed[i,k] = totFlow / rho[i,k]
						else:
							speed[i,k] = vMax[i]
						
						for j in range(len(id.paths)):
							if i in id.paths[j]:
								ttInst[j,k] += D / speed[i,k]
						
						if rho[i,k] > 0.1:
							headway[i,k] = ( (1/rho[i,k]) - L ) / speed[i,k] * 3600
						else:
							headway[i,k] = headway[i,k-1]
							
				print('TTS: ' + str(tts))
	return headway
				
				
# call the main function
try:
	mainOptimisation(True)
except mosek.Error as e:
	print("ERROR: %s" % str(e.errno))
	if e.msg is not None:
		print("\t%s" % e.msg)
		sys.exit(1)
except:
	import traceback
	traceback.print_exc()
	sys.exit(1)