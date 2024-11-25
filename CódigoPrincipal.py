#Bibliotecas utilizadas na implementação do código
import math
import time
import random
from collections import defaultdict

class MinHeap:
    def __init__(self):
        self.heap = []  #Lista para armazenar a ordem
        self.dict = {}  #Dicionário para armazenar os elementos

    def parent(self, i):
        # Determina o índice do pai de um vértice
        return (i - 1) // 2

    def left_child(self, i):
        # Determina o índice do filho esquerdo de um vértice
        return 2 * i + 1

    def right_child(self, i):
        # Determina o índice do filho direito de um vértice
        return 2 * i + 2

    def insert(self, key, value):
        # Adiciona um novo elemento ao heap
        self.dict[key] = value
        if key in self.heap:
            self.heap.remove(key)
        self.heap.append(key)
        self.heapify_up(len(self.heap) - 1)

    def heapify_up(self, i):
        # Correção da posição de um elemento no heap (para cima)
        while i != 0 and self.dict[self.heap[self.parent(i)]] > self.dict[self.heap[i]]:
            self.heap[self.parent(i)], self.heap[i] = self.heap[i], self.heap[self.parent(i)]
            i = self.parent(i)

    def extract_min(self):
        # Remove e retorna o menor elemento do heap
        if len(self.heap) == 0:
            return (None, None)
        
        # O menor valor é sempre a raiz (primeiro elemento)
        root_key = self.heap[0]
        
        # Coloca o último elemento na raiz e ajusta o heap
        last_key = self.heap.pop()
        if self.heap:
            self.heap[0] = last_key
            self.heapify_down(0)
        
        # Remove a entrada correspondente no dicionário
        min_value = self.dict.pop(root_key)
        
        return (root_key, min_value)

    def heapify_down(self, i):
        # Correção da posição de um elemento no heap (para baixo)
        smallest = i
        left = self.left_child(i)
        right = self.right_child(i)

        # Verifica se o filho esquerdo é menor que o nó atual
        if left < len(self.heap) and self.dict[self.heap[left]] < self.dict[self.heap[smallest]]:
            smallest = left
        
        # Verifica se o filho direito é menor que o nó atual
        if right < len(self.heap) and self.dict[self.heap[right]] < self.dict[self.heap[smallest]]:
            smallest = right
        
        # Se o menor não for o nó atual, faz a troca e continua
        if smallest != i:
            self.heap[i], self.heap[smallest] = self.heap[smallest], self.heap[i]
            self.heapify_down(smallest)

    def get_min(self):
        # Retorna o menor elemento sem removê-lo
        if len(self.heap) == 0:
            return None
        key = self.heap[0]
        return (key, self.dict[key])

    def size(self):
        # Determina o tamanho do heap
        return len(self.heap)

    def get_heap(self):
        # Retorna uma representação do heap
        return [(key, self.dict[key]) for key in self.heap]
            

class Grafo:
    def __init__(self, tipo, direc):
        #Características do grafo
        self.n=0 #Número de vértices                                                                     
        self.m=0 #Número de arestas
        self.tipo=tipo #Tipo de representação do grafo
        self.direc=direc #Grafo direcionado ou não
        self.arestas=[] #Lista de arestas do grafo  
        self.representacao=[] #Representação do grafo
        self.graus=[]  #Lista de graus dos vértices
        self.vertices=[] #Lista de vértices
        self.original = defaultdict(dict)  #Grafo original
        self.residual = defaultdict(dict)  #Grafo residual

        if self.tipo!='Matriz' and self.tipo!='Lista': #Erro para tipo de grafo inválido
            raise Exception('Tipo de grafo inválido!')
        if type(self.direc) != bool:  #Erro para entrada errada de direcionamento do grafo
            raise Exception('Entrada inválida!')
    
    def ImportarTxt(self, grafo_txt):
        #Função para ler um grafo de um arquivo texto
        with open(grafo_txt, 'r', encoding='utf-8') as arquivo:
            lista=[]
            for linha in arquivo:
                lista.append(linha.strip().split())
            self.n= int(lista[0][0])
            self.m= len(lista)-1
            self.arestas=lista[1:]
            self.graus=[0]*self.n
            if self.direc==False:   #Caso o grafo não seja direcionado
                for i in range(self.n):
                    if self.tipo=='Matriz':     #Criando a representação do grafo de acordo com o tipo escolhido pelo usuário
                        self.representacao.append([False]*self.n)
                    elif self.tipo=='Lista':    #Criando a representação do grafo de acordo com o tipo escolhido pelo usuário
                        self.representacao.append([])
                for aresta in self.arestas:
                    self.graus[int(aresta[0])-1]+=1
                    self.graus[int(aresta[1])-1]+=1
                    if self.tipo=='Matriz':
                        self.representacao[int(aresta[0])-1][int(aresta[1])-1]=int(aresta[2])
                        self.representacao[int(aresta[1])-1][int(aresta[0])-1]=int(aresta[2])
                        if int(aresta[2])<0:
                            raise Exception('A biblioteca ainda não implementa caminhos mínimos com pesos negativos.')
                    elif self.tipo=='Lista':
                        self.representacao[int(aresta[0])-1].append([int(aresta[1]), int(aresta[2])])
                        self.representacao[int(aresta[1])-1].append([int(aresta[0]), int(aresta[2])])
                        if int(aresta[2])<0:
                            raise Exception('A biblioteca ainda não implementa caminhos mínimos com pesos negativos.')
                if self.tipo=='Lista':
                    for i in range(self.n):
                        self.representacao[i]=sorted(self.representacao[i])

            if self.direc==True:   #Caso o grafo seja direcionado
                for i in range(self.n):
                    if self.tipo=='Matriz':     #Criando a representação do grafo de acordo com o tipo escolhido pelo usuário
                        self.representacao.append([False]*self.n)
                    elif self.tipo=='Lista':    #Criando a representação do grafo de acordo com o tipo escolhido pelo usuário
                        self.representacao.append([])
                for aresta in self.arestas:
                    self.graus[int(aresta[0])-1]+=1 #Salvando apenas o grau de saída de cada vértice
                    if self.tipo=='Matriz':
                        self.representacao[int(aresta[0])-1][int(aresta[1])-1]=int(aresta[2])
                        if int(aresta[2])<0:
                            raise Exception('A biblioteca ainda não implementa caminhos mínimos com pesos negativos.')
                    elif self.tipo=='Lista':
                        self.representacao[int(aresta[0])-1].append([int(aresta[1]), int(aresta[2])])
                        if int(aresta[2])<0:
                            raise Exception('A biblioteca ainda não implementa caminhos mínimos com pesos negativos.')
                    self.original[int(aresta[0])][int(aresta[1])]= {'capacity': int(aresta[2]), 'flow': 0}    #Criando o grafo original
                if self.tipo=='Lista':
                    for i in range(self.n):
                        self.representacao[i]=sorted(self.representacao[i])

    def residual_add_edge(self, u, v, capacity):
        #Função para adicionar uma aresta ao grafo
        self.residual[u][v] = {'capacity': capacity} 

    def create_residual_graph(self):
        #Função para criar um grafo residual
        for u in self.original:
            for v in self.original[u]:
                capacity = self.original[u][v]['capacity']
                flow = self.original[u][v]['flow']
                if capacity - flow > 0: #Aresta original
                    self.residual_add_edge(u, v, capacity - flow)
                if flow > 0: #Aresta reversa
                    self.residual_add_edge(v, u, flow)
        return self.residual 

    def BFS_Ford_Fulkerson(self, source, sink):
        #Função auxiliar para encontrar um caminho aumentante
        pai = [None] * self.n  
        nivel = [0] * self.n  
        marcados = [False] * self.n  
        fila = [source]  
        marcados[source-1] = True  

        while fila:
            v = fila.pop(0) 

            # Parar ao encontrar o vértice desejado
            if v == sink:
                path = []
                current = sink
                while current != source:
                    prev = pai[current - 1]
                    path.append((prev, current))
                    current = prev
                path.reverse()
                return path

            for w in self.residual[v]:
                if not marcados[w - 1] and self.residual[v][w]['capacity'] > 0:
                    fila.append(w)
                    marcados[w - 1] = True
                    pai[w - 1] = v
                    nivel[w - 1] = nivel[v - 1] + 1

        # Nenhum caminho aumentante encontrado
        return None 

    def update_flow(self, path, bottleneck):
        #Função auxiliar para atualizar os grafos original e residual com o novo gargalo

        for u, v in path:
        # Atualizar fluxo no grafo original
            if v in self.original[u]:  
                self.original[u][v]['flow'] += bottleneck

        # Atualizar valores no grafo residual
            if v in self.residual[u]:
                self.residual[u][v]['capacity'] -= bottleneck
                if self.residual[u][v]['capacity'] == 0:
                    del self.residual[u][v]  # Desconsiderar arestas sem capacidade

            if u not in self.residual[v]:
                self.residual[v][u] = {'capacity': bottleneck}
            else:
                self.residual[v][u]['capacity'] += bottleneck

    def calculate_bottleneck(self, path):
        #Função auxiliar para calcular o gargalo do caminho
        return min(self.residual[u][v]['capacity'] for u, v in path)
    
    def Ford_Fulkerson(self, source, sink):
        #Função para determinar o fluxo máximo de um grafo, por meio do algoritmo de Ford-Fulkerson
        max_flow = 0
        self.create_residual_graph()

        while True:
            path = self.BFS_Ford_Fulkerson(source, sink)
            if not path:
                break  # Não há mais caminhos aumentantes

            bottleneck = self.calculate_bottleneck(path)
            self.update_flow(path, bottleneck)
            max_flow += bottleneck

        with open('InfoFordFulkerson.txt', 'w', encoding='utf-8') as infoFordFulkerson:
            infoFordFulkerson.write(f'Fluxo Máximo: {max_flow}\n')
            infoFordFulkerson.write('\n')
            infoFordFulkerson.write('Fluxo nas Arestas:\n')
            infoFordFulkerson.write('\n')
            for u in self.original:
                for v in self.original[u]:
                    infoFordFulkerson.write(f'Aresta ({u}, {v}): {self.original[u][v]["flow"]}\n')

        return max_flow, self.original

    def BFS(self, raiz):
        #Função para realizar uma busca em largura no grafo
        pai=[None]*self.n
        nivel=[0]*self.n
        marcados=[False]*self.n
        fila=[raiz]
        marcados[raiz-1]=True
        while len(fila)!=0:
            v=fila.pop(0)
            if self.tipo=='Lista':      #Adaptação ao tipo escolhido pelo usuário
                for w in self.representacao[v-1]:
                    if not marcados[w[0]-1]:
                        fila.append(w[0])
                        marcados[w[0]-1]=True
                        pai[w[0]-1]=v
                        nivel[w[0]-1]=nivel[v-1]+1
            if self.tipo=='Matriz':     #Adaptação ao tipo escolhido pelo usuário
                for w in range(self.n):
                    if type(self.representacao[v-1][w])!=bool:
                        if not marcados[w]:
                            fila.append(w+1)
                            marcados[w]=True
                            pai[w]=v
                            nivel[w]=nivel[v-1]+1
        with open('InfoBFS.txt', 'w', encoding='utf-8') as infoBFS: #Criação de um arquivo com as informações da busca em largura
            infoBFS.write('Informações da Busca em Largura:\n')
            infoBFS.write('\n')
            for i in range(self.n):
                infoBFS.write(f'Vértice {i+1}: Pai: {pai[i]}, Nível: {nivel[i]}\n')

    def DFS(self,raiz):
        #Função para realizar uma busca em profundidade no grafo
        pai=[None]*self.n
        nivel=[0]*self.n
        marcados=[False]*self.n
        pilha=[raiz]
        while len(pilha)!=0:
            v=pilha.pop()
            if not marcados[v-1]:
                marcados[v-1]=True
                if self.tipo=='Lista':      #Adaptação ao tipo escolhido pelo usuário
                    for w in sorted(self.representacao[v-1], reverse=True):
                        pilha.append(w[0])
                        if not marcados[w[0]-1]:
                            pai[w[0]-1]=v
                            nivel[w[0]-1]=nivel[v-1]+1
                elif self.tipo=='Matriz':       #Adaptação ao tipo escolhido pelo usuário
                    for w in range(self.n-1, -1, -1):
                        if type(self.representacao[v-1][w])!=bool:
                            pilha.append(w+1)
                            if not marcados[w]:
                                pai[w]=v
                                nivel[w]=nivel[v-1]+1
        with open('InfoDFS.txt', 'w', encoding='utf-8') as infoDFS: #Criação de um arquivo com as informações da busca em profundidade
            infoDFS.write('Informações da Busca em Profundidade:\n')
            infoDFS.write('\n')
            for i in range(self.n):
                infoDFS.write(f'Vértice {i+1}: Pai: {pai[i]}, Nível: {nivel[i]}\n')

    def BFS_Mais_Distante(self, raiz):
        #Função para realizar uma busca em largura no grafo e determinar o vértice mais distante da raiz
        #Função elaborada para auxiliar outras funções
        pai=[None]*self.n
        Mais_Distante=raiz
        nivel=[0]*self.n
        marcados=[False]*self.n
        fila=[raiz]
        marcados[raiz-1]=True
        while len(fila)!=0:
            v=fila.pop(0)
            if self.tipo=='Lista':
                for w in sorted(self.representacao[v-1]):
                    if not marcados[w[0]-1]:
                        fila.append(w[0])
                        marcados[w[0]-1]=True
                        pai[w[0]-1]=v
                        nivel[w[0]-1]=nivel[v-1]+1
                        if nivel[w[0]-1]>nivel[Mais_Distante-1]:
                            Mais_Distante=w[0]
            if self.tipo=='Matriz':
                for w in range(self.n):
                    if type(self.representacao[v-1][w])!=bool:
                        if not marcados[w]:
                            fila.append(w+1)
                            marcados[w]=True
                            pai[w]=v
                            nivel[w]=nivel[v-1]+1
                            if nivel[w]>nivel[Mais_Distante-1]:
                                Mais_Distante=w
        return pai, nivel, Mais_Distante
    
    def distancia(self, x, y):
        #Função para calcular a distância entre dois vértices
        distancia=self.Djikstra_Return(x)
        if distancia[y-1]==math.inf:
            return 'Não existe caminho entre os vértices!' 
        else:
            return distancia[y-1]

    def diametro(self):
        #Função para calcular o diâmetro do grafo
        lista_auxiliar=[]
        for i in range(1, self.n+1):
            distancia=self.Djikstra_Return(i)
            lista_auxiliar.append(max(distancia))
        diametro= max(lista_auxiliar)
        return diametro    

    def ComponentesConexas(self):
        #Função para determinar as componentes conexas do grafo
        ComponentesConexas=[]
        Marcados=[]
        Num_Componentes=0
        pai, nivel, Mais_Distante= self.BFS_Mais_Distante(1)
        for i in range(self.n):
            if len(Marcados)==self.n:
                break
            elif pai[i]== None and i+1 not in Marcados:
                Num_Componentes+=1
                ComponentesConexas.append([i+1])
                Marcados.append(i+1)
                pai, nivel, Mais_Distante= self.BFS_Mais_Distante(i+1)
                for j in range(len(pai)):
                    if pai[j]!=None:
                        Marcados.append(j+1)
                        ComponentesConexas[Num_Componentes-1].append(j+1)
        resposta= f'Número de Componentes Conexas: {Num_Componentes} \n' 
        for z in range(Num_Componentes): 
            resposta+= f'Componente {z+1}: Tamanho: {len(sorted(ComponentesConexas, key=len, reverse=True)[z])}, Vértices: {sorted(ComponentesConexas, key=len, reverse=True)[z]}\n'
        return resposta
    
    def DiametroAprox(self):
        #Função para calcular um diâmetro aproximado do grafo
        vertice=random.randint(1, self.n)
        distancia=self.Djikstra_Return(vertice)      #Encontrando um vértice afastado do grafo
        distancia2=self.Djikstra_Return(distancia.index(max(distancia)))    #Aproximando o diâmetro a partir do vértice mais distante de um vértice afastado do grafo
        DiametroAprox= max(distancia2)
        return DiametroAprox
    
    def Djikstra(self, raiz):
        #Função para determinar o menor caminho entre um vértice e todos os outros do grafo, implementada com o uso de um vetor
        distancias=[math.inf]*self.n
        marcados=[False]*self.n
        distancias[raiz-1]=0
        caminhos=[]
        for i in range(self.n):
            caminhos.append([])
        while False in marcados:
            menor=math.inf
            for j in range(self.n):
                if not marcados[j] and distancias[j]<menor:
                    menor=distancias[j]
                    vertice=j+1
            marcados[vertice-1]=True
            if self.tipo=='Lista':
                for w in self.representacao[vertice-1]:
                    if distancias[vertice-1]+w[1]<distancias[w[0]-1]:
                        distancias[w[0]-1]=distancias[vertice-1]+w[1]
                        caminhos[w[0]-1]= caminhos[vertice-1] + [vertice]
            if self.tipo=='Matriz':
                for w in range(self.n):
                    if type(self.representacao[vertice-1][w])!=bool:
                        if distancias[vertice-1]+self.representacao[vertice-1][w]<distancias[w]:
                            distancias[w]=distancias[vertice-1]+self.representacao[vertice-1][w]
                            caminhos[w]= caminhos[vertice-1] + [vertice]
        with open('InfoDjikstra.txt', 'w', encoding='utf-8') as infoDjikstra: #Criação de um arquivo com as informações do menor caminho
            infoDjikstra.write('Menor Caminho:\n')
            infoDjikstra.write('\n')
            for i in range(self.n):
                infoDjikstra.write(f'Vertice {raiz} para Vértice {i+1}: {distancias[i]}\n')
                infoDjikstra.write(f'Caminho: {caminhos[i]} -> {i+1}\n')

    def DjikstraHeap(self, raiz):
        #Função para determinar o menor caminho entre um vértice e todos os outros do grafo, implementada com auxílio de um Heap de mínimos
        distancias=[math.inf]*self.n
        distanciasHeap=MinHeap()
        marcados=[False]*self.n
        distancias[raiz-1]=0
        distanciasHeap.insert(raiz, 0)
        caminhos=[]
        for i in range(self.n):
            caminhos.append([])
        while distanciasHeap.size()!=0:
            vertice, valor= distanciasHeap.extract_min()
            marcados[vertice-1]=True
            if self.tipo=='Lista':
                for w in self.representacao[vertice-1]:
                    if distancias[vertice-1]+w[1]<distancias[w[0]-1]:
                        distancias[w[0]-1]=distancias[vertice-1]+w[1]
                        distanciasHeap.insert(w[0], distancias[vertice-1]+w[1])
                        caminhos[w[0]-1]= caminhos[vertice-1] + [vertice]
            if self.tipo=='Matriz':
                for w in range(self.n):
                    if type(self.representacao[vertice-1][w])!=bool:
                        if distancias[vertice-1]+self.representacao[vertice-1][w]<distancias[w]:
                            distancias[w]=distancias[vertice-1]+self.representacao[vertice-1][w]
                            distanciasHeap.insert(w+1, distancias[vertice-1]+self.representacao[vertice-1][w])
                            caminhos[w]= caminhos[vertice-1] + [vertice]
        with open('InfoDjikstraHeap.txt', 'w', encoding='utf-8') as infoDjikstraHeap: #Criação de um arquivo com as informações do menor caminho
            infoDjikstraHeap.write('Menor Caminho:\n')
            infoDjikstraHeap.write('\n')
            for i in range(self.n):
                infoDjikstraHeap.write(f'Vertice {raiz} para Vértice {i+1}: {distancias[i]}\n')
                infoDjikstraHeap.write(f'Caminho: {caminhos[i]} -> {i+1}\n')

    def Djikstra_Return(self, raiz):
        #Função para determinar o menor caminho entre um vértice e todos os outros do grafo
        #Função auxiliar de outras funções
        distancias=[math.inf]*self.n
        distanciasHeap=MinHeap()
        marcados=[False]*self.n
        distancias[raiz-1]=0
        distanciasHeap.insert(raiz, 0)
        while distanciasHeap.size()!=0:
            vertice, valor= distanciasHeap.extract_min()
            marcados[vertice-1]=True
            if self.tipo=='Lista':
                for w in self.representacao[vertice-1]:
                    if distancias[vertice-1]+w[1]<distancias[w[0]-1]:
                        distancias[w[0]-1]=distancias[vertice-1]+w[1]
                        distanciasHeap.insert(w[0], distancias[vertice-1]+w[1])
            if self.tipo=='Matriz':
                for w in range(self.n):
                    if type(self.representacao[vertice-1][w])!=bool:
                        if distancias[vertice-1]+self.representacao[vertice-1][w]<distancias[w]:
                            distancias[w]=distancias[vertice-1]+self.representacao[vertice-1][w]
                            distanciasHeap.insert(w+1, distancias[vertice-1]+self.representacao[vertice-1][w])
        return distancias

    def CriarTxt(self):
        #Função para criar um arquivo texto de saída com informações do grafo
        with open('Info.txt', 'w', encoding='utf-8') as info:
            if self.n%2==0:
                grau_med= (sorted(self.graus)[(self.n//2)-1]+sorted(self.graus)[(self.n//2)])/2
            else:
                grau_med= sorted(self.graus)[(math.ceil(self.n/2))-1]
            info.write('Informações do Grafo:\n')
            info.write('\n')
            info.write(f'Número de Vértices: {self.n}\n')
            info.write(f'Número de Arestas: {self.m}\n')
            info.write(f'Grau Mínimo: {min(self.graus)}\n')
            info.write(f'Grau Máximo: {max(self.graus)}\n')
            info.write(f'Grau Médio: {sum(self.graus)/self.n}\n')
            info.write(f'Grau Mediano: {grau_med}\n')
            info.write('\n')
            info.write('Informações das Componentes Conexas:\n')
            info.write(self.ComponentesConexas())	

