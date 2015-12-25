#coding=utf-8
'''
实现隐马尔科夫模型的基本方法，
输入为状态转移矩阵，观测矩阵，初始状态概率向量，观测序列
实现前向算法和后向算法计算观测序列出现的概率
实现维特比算法找当前观测序列下最可能的状态序列
实现在给定模型和观测下，t时刻处于p状态的概率，t,p在main函数中指定
'''
from numpy import *

class HMM:

	def __init__(self):
		self.A=array([(0.5,0.2,0.3),(0.3,0.5,0.2),(0.2,0.3,0.5)])
		self.B=array([(0.5,0.5),(0.4,0.6),(0.7,0.3)])
		self.pi=array([(0.2),(0.4),(0.4)])
		self.o=[0,1,0]
		self.t=len(self.o)#观测序列长度
		self.m=len(self.A)#状态集合个数
		self.n=len(self.B[0])#观测集合个数

	def qianxiang(self):
		#t时刻部分观测序列为o1,o2,ot,状态为qi的概率用矩阵x表示，
		#则矩阵大小行数为观测序列数，列数为状态个数
		self.x=array(zeros((self.t,self.m)))
		#先计算出时刻1时，观测为o1,状态为qi的概率
		for i in range(self.m):
			self.x[0][i]=self.pi[i]*self.B[i][self.o[0]]
		for j in range(1,self.t):
			for i in range(self.m):
				#前一时刻所有状态的概率乘以转移概率得到i状态概率
				#i状态的概率*i状态到j观测的概率
				temp=0
				for k in range(self.m):
					temp=temp+self.x[j-1][k]*self.A[k][i]
				self.x[j][i]=temp*self.B[i][self.o[j]]
		result=0
		for i in range(self.m):
			result=result+self.x[self.t-1][i]
		print u"前向概率矩阵及当前观测序列概率如下："
		print self.x
		print result

	def houxiang(self):
		#t时刻状态为qi,从t+1到T观测为ot+1,ot+2,oT的概率用矩阵y表示
		#则矩阵大小行数为观测序列数，列数为状态个数
		self.y=array(zeros((self.t,self.m)))
		#下面为对最终时刻的所有状态，接下来的观测序列概率初始化为1
		#(可以理解为接下来没有观测所有为1)
		for i in range(self.m):
			self.y[self.t-1][i]=1
		j=self.t-2
		#j时刻为i状态，转移到k状态，k状态观测为oj+1,
		#再乘以j+1时刻状态为k的B矩阵的值，对k遍历相加，
		#即为j时刻i状态前提下，后面满足观测序列的概率
		while(j>=0):
			for i in range(self.m):
				for k in range(self.m):
					self.y[j][i]+=self.A[i][k]*self.B[k][self.o[j+1]]*self.y[j+1][k]
			j=j-1
		#第一个状态任意，观测为o1,所以对所有第一个状态概率相加
		result=0
		for i in range(self.m):
			result=result+self.pi[i]*self.B[i][self.o[0]]*self.y[0][i]
		print u'后向概率矩阵及当前观测序列概率如下：'
		print self.y
		print result

	def get_stateprobability(self,t,p):
		#打印在观测为self.o的前提下，t时刻，处于状态p的概率,
		#self.x[t][p]表示到t时刻观测为o1,o2,ot,状态为p的概率
		#self.y[t][p]表示在t时刻状态为p的前提下，接下来观测为ot+1,ot+2,oT的概率
		#self.x[t][p]*self.y[t][p]即表示观测为self.o，且t时刻处于状态p的概率,
		#利用贝叶斯公式，除以观测为self.o的概率即为所求
		if(t>self.t or p>self.m):
			print u'输入数据超过范围'
			return
		print u'在时刻'+str(t)+u'处于状态'+str(p)+u'的概率是：'
		temp=self.x[t-1][p-1]*self.y[t-1][p-1]
		total=0
		for i in range(self.m):
			total=total+self.x[t-1][i]*self.y[t-1][i]
		print temp/total

	def viterbi(self):
		#利用模型和观测序列找出最优的状态序列
		#时刻t时，很多路径可以到达状态i,且观测为self.o,
		#每个路径都有自己的概率，最大的概率用矩阵z记录,前一个状态用d矩阵记录
		self.z=array(zeros((self.t,self.m)))
		self.d=array(zeros((self.t,self.m)))
		for i in range(self.m):
			self.z[0][i]=self.pi[i]*self.B[i][self.o[0]]
			self.d[0][i]=0
		for j in range(1,self.t):
			for i in range(self.m):
				maxnum=self.z[j-1][0]*self.A[0][i]
				node=1
				for k in range(1,self.m):
					temp=self.z[j-1][k]*self.A[k][i]
					if(maxnum<temp):
						maxnum=temp
						node=k+1
				self.z[j][i]=maxnum*self.B[i][self.o[j]]
				self.d[j][i]=node
		#找到T时刻概率最大的路径
		max_probability=self.z[self.t-1][0]
		last_node=[1]
		temp=0
		for i in range(1,self.m):
			if(max_probability<self.z[self.t-1][i]):
				max_probability=self.z[self.t-1][i]
				last_node[0]=i+1
				temp=i
		i=self.t-1
		#self.d[t][p],t时刻状态为p的时候，t-1时刻的状态
		while(i>=1):
			last_node.append(self.d[i][temp])
			i=i-1
		temp=['o']
		temp[0]=int(last_node[len(last_node)-1])
		j=len(last_node)-2
		while j>=0:
			temp.append(int(last_node[j]))
			j=j-1
		print u'路径节点分别为'
		print temp
		print u'该路径概率为'+str(max_probability)


if __name__ == '__main__':
	test=HMM()
	test.qianxiang()
	test.houxiang()
	test.get_stateprobability(3,3)
	test.viterbi()
