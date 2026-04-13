# TODO: 1) понять какой этап при вызове занимает так долго
# 2) как сохранять переменные для символьных функий? -> ссылки на экземпляры классов -> 


from Lego import *
import math as m
from ModelFiles import to_logOut
import inspect
import Lego_smbFun as smb

def simple_split(text):
    parts = []
    part = ""
    depth = 0
    
    for ch in text:
        if ch == '(':
            depth += 1
            part += ch
        elif ch == ')':
            depth -= 1
            part += ch
        elif ch == ',' and depth == 0:
            parts.append(part.strip())
            part = ""
        else:
            part += ch
    
    if part:
        parts.append(part.strip())
    
    return parts

def parse_func_args(text):
    """
    Парсит строку вида: 
    - "Spwl(6, 10)" -> ("Spwl", ["6", "10"])
    - "Spwl()" -> ("Spwl", [])
    - "Pwl(1, 5)" -> ("Pwl", ["1", "5"])
    """
    text = text.strip()
    
    # Ищем первую открывающую скобку
    name_end = text.find('(')
    
    # Если скобок нет
    if name_end == -1:
        return text, []
    
    # Имя функции - часть до скобки
    func_name = text[:name_end].strip()
    
    # Часть внутри скобок (убираем последнюю закрывающую скобку)
    args_text = text[name_end+1:-1].strip() if text.endswith(')') else text[name_end+1:]
    
    # Разбираем аргументы с учетом вложенности
    args = args_text.split(",")
    
    return func_name, args


class MixedFun (smb.smbFun) :
    """
    Синтаксис: f ( X, V ) = [SPWL,Mesh,Polynome(6)];
    """

    def __init__(self, Vname='', As=[], discr_str=None, param=False, Finitialize=0, DataReadFrom='', Data=[], Type='', Domain=None, ReadFrom=''):
        smb.smbFun.__init__(self, Vname, As, param, -1, Finitialize, DataReadFrom, Data, Type, Domain, ReadFrom)
        self.SymbolInteg=False
        self.SymbolDiffer =False
        self.Deriv1=False
        self.do_print = True
        self.do_print_ext = False
        self.ArgNorm = True

        self.discr_str = discr_str
        self.discr_list = dict()
        self.smb_discr_list = dict()
        self.grd_discr_list = dict()
        self.include_smb = True
        self.include_grd = True
        self.pre_smbF = None
        self.smbF = self.complex_smbF
        self.DERIV2 = self.DERIV2_mine
        self.F = self.F_mine
        self.Fijk = self.Fijk_mine
        self.var_to_grd = self.var_to_grd_mine
        self.grd_norm = "Node" # "Node" or "Norm01"
        self.smb_norm = "Norm01" # "Node" or "Norm01"

        self.process_discr_string()
        self.smb_axis = [int(_) for _ in self.smb_discr_list.keys()]
        self.grd_axis = [int(_) for _ in self.grd_discr_list.keys()]
        self.print_info()
        # self.discr_types = discr_types
        self.dim = len(self.discr_list)

        self.discr_methods = { 
            'SPWL': self._SPWL_discr,
            'Mesh': self._Mesh_discr}
        
        # self.discr_methods = { 
        #     'SPWL': self._SPWL_discr,
        #     'Polynome': self._polynomial_discr,
        #     'Fourier': self._fourier_discr,
        #     'SOS2': self._SOS2_discr}
    
    def process_discr_string(self):
        # пусть на вход входит "[SPWL,Mesh,Polynome(6)]"
        discr_list = simple_split(self.discr_str[1:-1])
        for i, discr_method in enumerate(discr_list):
            if discr_method.startswith('SPWL') == 1 or discr_method.startswith('Mesh') == 1:
                self.grd_discr_list[str(i)] = discr_method
            else:
                self.smb_discr_list[str(i)] = discr_method

        self.discr_list = discr_list
        if len(self.smb_discr_list) == 0: self.include_smb = False
        if len(self.grd_discr_list) == 0: self.include_grd = False

        return 0
    
    def print_info(self):
        print("Function configuration:")
        print(f"    grid types = {self.grd_discr_list}, grid_coords = {self.grd_axis}")
        print(f"    symbolic types = {self.smb_discr_list}, smb_coords = {self.smb_axis}")

    
    @staticmethod
    def _CHKS_smoothing(z, e): # CHKS smoothing
        return 0.5 * (z + py.sqrt(z**2 + e))
    

    def _SPWL_discr(self, x_points, y_points, x_eval, epsilon=SvF.Epsilon): 
        """
        SPWL-дискретизация, Pyomo-only версия.

        x_points : Set object
            Множество узлов по оси
        y_points : IndexedVar / IndexedExpression
            Значения функции в узлах
        x_eval : Pyomo Expression
            Точка, в которой считаем интерполяцию
        epsilon : float
            Параметр сглаживания
        """
        #print(x_points)
        if self.do_print: print(f"Внутри SPWL_discr с x_eval = {x_eval}\n")       
        epsilon = float(epsilon)
        step = 1
        A_prev = 0
        A_curr = (y_points[1] - y_points[0]) / step
        result = y_points[0] + (A_curr - A_prev) * (x_eval - x_points[0])

        for i in range(2, len(x_points)):
            A_next = (y_points[i] - y_points[i - 1]) / step

            result += (A_next - A_curr) * self._CHKS_smoothing(
                x_eval - x_points[i - 1], epsilon)

            A_curr = copy(A_next)

        return result
    

    def _SPWL_discr_old(self, x_points, y_points, x_eval, epsilon=SvF.Epsilon): # без приведения к узлам вида 0, 1, 2....
        """
        SPWL-дискретизация, Pyomo-only версия.

        x_points : Set object
            Множество узлов по оси
        y_points : IndexedVar / IndexedExpression
            Значения функции в узлах
        x_eval : Pyomo Expression
            Точка, в которой считаем интерполяцию
        epsilon : float
            Параметр сглаживания
        """
        epsilon = float(epsilon)
        step = x_points.step
        A_prev = 0

        A_curr = (y_points[1] - y_points[0]) / step
        result = y_points[0] + (A_curr - A_prev) * (x_eval - x_points.Val[0])

        for i in range(2, x_points.Ub+1):
            A_next = (y_points[i] - y_points[i - 1]) / step

            result += (A_next - A_curr) * self._CHKS_smoothing(
                x_eval - x_points.Val[i - 1], epsilon)

            A_curr = copy(A_next)

        return result
    

    @staticmethod
    def _Mesh_discr(x_points, y_points, x_eval): # для сетки, нормированной по узлам, Node
        # Находим интервал
        x_i = x_points[0]
        x_i1 = x_points[1]

        # print(f"В mesh x_eval = {x_eval}")
        if x_i <= 0: 
            return y_points[0]
    
        if x_i >= len(x_points)-1: 
            return y_points[len(x_points)-1]
        
        # if x_i == len(x_points)-2: x_i = x_i - 1
        #print(f"x_i = {x_i}", f"x_i1 = {x_i1}")
        # Вычисляем веса для двух соседних точек 
        step = 1
        lambda_i = (x_i1 - x_eval) / step  # вес для i-й точки
        lambda_i1 = (x_eval - x_i) / step  # вес для (i+1)-й точки
        
        # Значение интерполяции
        value = lambda_i * y_points[x_i] + lambda_i1 * y_points[x_i+1]
        
        return value 
    

    # @staticmethod
    # def _Mesh_discr_old(x_points, y_points, x_eval): # без приведения к узлам вида 0, 1, 2....
    #     # Находим интервал
    #     i = x_points.ValToInd(x_eval)
        
    #     # Точное попадание в узел
    #     # if np.isclose(x_eval, x_points[i]):
    #     #     return y_points[i] #, [i], [1.0]
        
    #     if i <= 0: 
    #         return y_points[0]
    
    #     if i >= x_points.Ub: 
    #         return y_points[x_points.Ub]
        
    #     # Вычисляем веса для двух соседних точек (SOS2)
    #     x_i, x_i1 = x_points.Val[i], x_points.Val[i+1]
    #     d_total = x_points.step
        
    #     lambda_i = (x_i1 - x_eval) / d_total  # вес для i-й точки
    #     lambda_i1 = (x_eval - x_i) / d_total  # вес для (i+1)-й точки
        
    #     # Значение интерполяции
    #     value = lambda_i * y_points[i] + lambda_i1 * y_points[i+1]
    #     # to_logOut("Bruh")
        
    #     return value
    

    def interpoleNode ( self, argNode, lev=None, to_Node=False):
        if lev is None :
            lev = len(self.grd_discr_list)

            # приведение сеточной части координат к узлам
            if to_Node == True: 
                argNode1 = self.real_to_normalized(copy(argNode), mode_grd="Node", mode_smb="Real")
            else:
                argNode1 = copy(argNode)

        else:
            argNode1 = copy(argNode)

        if lev > 0 :
            method_name = list(self.grd_discr_list.values())[lev-1] #[0]
            axis_num = int(list(self.grd_discr_list.keys())[lev-1])
            # print(f"axis_num = {axis_num}, размер - {self.A[axis_num].NodS}") # debug
            method_name, extra_args = parse_func_args(method_name)
            discr_func = self.discr_methods.get(method_name)

            if discr_func is None:
                raise ValueError(
                    f"Неизвестный метод дискретизации '{method_name}'. "
                    f"Доступные: {list(self.discr_methods.keys())}"
            )
            #print(f"текущий набор индексов {indices} и уровень {lev}") # debug
            # print(f"Набор узлов: {self.A[0].Val}, {self.A[1].Val}") # debug
            if self.do_print_ext:
                print(f"На этапе {lev} вызывается {discr_func.__name__} c аргументами {extra_args}", "argNode:") # debug
                [print(argNode1[i]) for i in range(len(argNode1))]

            if method_name == "Mesh":
                index = argNode[axis_num]
                args_list = [[self.A[axis_num].Val[int(floor(index))], self.A[axis_num].Val[int(ceil(index))]], 
                            [self.interpoleNode(argNode1[:axis_num] + [i] + argNode1[axis_num+1:], lev-1) for i in [int(floor(index)), int(ceil(index))]], 
                            argNode1[axis_num]]
            else:
                args_list = [list(self.A[axis_num].NodS), 
                            [self.interpoleNode(argNode1[:axis_num] + [i] + argNode1[axis_num+1:], lev-1) for i in self.A[axis_num].NodS], 
                            argNode1[axis_num]]
            # args_list = [self.A[axis_num], [self.interpolNode([i] + indices, lev-1) 
            #                                for i in self.A[axis_num].NodS], argNode1[axis_num]]
            if extra_args: args_list.extend(extra_args)

            if self.do_print_ext:
                print(f"перед вызовом функции интерполяции, argNode={argNode}, lev={lev}, axis_num={axis_num}") # debug

            ret = discr_func(*args_list)
            
            # print(f"self.A[lev-1] = {self.A[lev-1].Val}")
            # print(f"[self.interpolNode([i] + indices, lev-1) \
            #                                for i in self.A[lev-1].NodS] = {[self.interpolNode([i] + indices, lev-1) \
            #                                for i in self.A[lev-1].NodS]}")
            
            return ret
        
        else :
            if self.do_print_ext:
                print(f"В interpole на этапе {lev} подставляется значение", "argNode:", end=' ') # debug
                print("[", end='')
                [print(argNode1[i], end=',') for i in range(len(argNode1))]
                print("]\n", end='')

            return self.fetch_value(argNode1)


    def fetch_value(self, argNode):
        if self.include_smb == False:
            return self.gr[tuple(argNode)]
        
        else:
            return self.smbF(argNode, "Don't denorm it")
            # if grd_axis:
            #     self.smbF(smb_cords, self.gr[tuple(grd_cords)])
            # else:
            #     self.smbF(smb_cords, [])


    def distribute_args(self, argNode):
        smb_coords = [argNode[_] for _ in self.smb_axis]
        grd_coords = [argNode[_] for _ in self.grd_axis]

        # print(f"smb_axis = {self.smb_axis}")
        # print(f"gr_axis = {self.grd_axis}")
        #if self.do_print: print(f"grd_coords = {grd_coords}, smb_coords = {smb_coords}")
        return grd_coords, smb_coords
         

    def complex_smbF(self, argNode, flag="do"):
        if self.include_smb == False:
            caller_name = inspect.currentframe().f_back.f_code.co_name
            if self.do_print: print(f"E_pre_smbF00 вызвана из: {caller_name}()")
            #argNode = self.normalized_grd_to_node(argNode, flag)
            if self.do_print: print(f"argNode = {argNode}, type = {type(argNode)}")
            return self.gr[tuple(argNode)]
        else:
            grd_coords, smb_coords = self.distribute_args(argNode)
            caller_name = inspect.currentframe().f_back.f_code.co_name
            if self.do_print: print(f"E_pre_smbF00 вызвана из: {caller_name}()")
            if self.include_grd == True:
                pass
                #grd_coords = self.normalized_grd_to_node(grd_coords, flag)

            return self.pre_smbF(smb_coords, tuple(grd_coords))


    def Ftbl ( self, n ) : 
        caller_name = inspect.currentframe().f_back.f_code.co_name
        if self.do_print: print(f"E_pre_smbF00 вызвана из: {caller_name}()")
        argNode = [a.dat[n] for a in self.A]
        if self.do_print: print('In Ftbl: n and argNode: ', n, argNode)
        argNode = self.real_to_normalized(argNode)
        if self.do_print: print('In Ftbl: n and argNode: ', n, argNode)
        if self.include_grd == False: return self.smbF(argNode)
        else: return self.interpoleNode(argNode)
    
#     def F ( self, ArS_real ) :  #   real args, по идее ненужный метод
#         print("F evoked "*5) # не бывает такого, видимо 0_0
#    #     if self.ArgNorm:  SvF.F_Arg_Type = 'N'     #  заплатка для ArgNorm для fNi_fon(X,Y) символ функции   Ni(X,Y)  = Ni_fon(X,Y) + fon
#         ret = self.smbF(self.real_to_normalized(ArS_real)) + self.V.avr    # gr = self.grd
#     #    SvF.F_Arg_Type = ''
#         return ret
    

    def real_to_normalized ( self, ArS_real, mode_grd=None, mode_smb=None) : 
        new_ArS = list() 
        if mode_grd == None: mode_grd = self.grd_norm
        if mode_smb == None: mode_smb = self.smb_norm

        for i, el in enumerate(ArS_real):
            if i in self.grd_axis:
                if mode_grd == "Norm01":
                    new_ArS.append((el-self.A[i].min)/self.A[i].ma_mi)
                elif mode_grd == "Node":
                    new_ArS.append((el-self.A[i].min)/self.A[i].step) 
                elif mode_grd == "Real":
                    new_ArS.append(el)
                else:
                    print("[WARNING] No normalization for grd")
                    new_ArS.append(el)
            else:
                if mode_smb == "Norm01":
                    new_ArS.append((el-self.A[i].min)/self.A[i].ma_mi) 
                elif mode_smb == "Node":
                    new_ArS.append((el-self.A[i].min)/self.A[i].step)   
                elif mode_smb == "Real":
                    new_ArS.append(el)
                else:
                    print("[WARNING] No normalization for smb")
                    new_ArS.append(el)     

        return new_ArS
    

    def node_to_normalized ( self, ArS_real, mode_grd=None, mode_smb=None) : 
        new_ArS = list() 
        if mode_grd == None: mode_grd = self.grd_norm
        if mode_smb == None: mode_smb = self.smb_norm

        for i, el in enumerate(ArS_real):
            if i in self.grd_axis:
                if mode_grd == "Norm01":
                    new_ArS.append(el*self.A[i].step/self.A[i].ma_mi)
                elif mode_grd == "Node":
                    pass
                elif mode_grd == "Real":
                    new_ArS.append((el + self.A[i].min)*self.A[i].step)
                else:
                    print("[WARNING] No normalization for grd")
                    new_ArS.append(el)
            else:
                if mode_smb == "Norm01":
                    new_ArS.append(el*self.A[i].step/self.A[i].ma_mi)
                elif mode_smb == "Node":
                    pass
                elif mode_smb == "Real":
                    new_ArS.append((el + self.A[i].min)*self.A[i].step)
                else:
                    print("[WARNING] No normalization for smb")
                    new_ArS.append(el)     

        return new_ArS
    

    def normalized_grd_to_node ( self, args, flag) : 
        new_args = list(copy(args))
        if self.do_print: print("was, ", new_args)
        #if self.grd_norm == "Norm01" and flag == "do":
        if flag == "do":
            for i, el in enumerate(args):
                if self.do_print: print("Внутри обратной нормировки, step = ", self.A[self.grd_axis[i]].step)
                new_args[i] = el / self.A[self.grd_axis[i]].step * self.A[self.grd_axis[i]].ma_mi
                #print("step, ", self.A[self.grd_axis[i]].step)
                #if abs(round(new_args[i]) - new_args[i]) < 0.01: new_args[i] = round(new_args[i]) ### IMPORTANT IMPORTANTIMPORTANTIMPORTANTIMPORTANT!!!!!!!!!!

        # new_args = [int(i) for i in new_args]   ### IMPORTANT IMPORTANTIMPORTANTIMPORTANTIMPORTANTIMPORTANT!!!!!!!!!!
        if self.do_print: print("became, ", new_args)
        return new_args

        # if self.ArgNorm  :  return  [(ArS_real[i]-a.min)/a.ma_mi for i, a in enumerate(self.A)]
        # else                    :  return ArS_real


    def Fijk_mine ( self, ijk ) :  # считает значение в узле с данными из таблицы по индексам
        arg_real = [ a.Val[ijk[ia]] for ia, a in enumerate (self.A) ]
        if self.do_print: print(f"arg_real = {arg_real}, arg_normalized = {self.real_to_normalized(arg_real)}")
        ret = self.interpoleNode(self.real_to_normalized(arg_real))
        if self.do_print: print(ret + self.V.avr) 
        return ret + self.V.avr    # gr = self.grd


    def var_to_grd_mine (self) :
            print("--- var_to_grd evoked")
            if self.dim == 1:
                for i in range(self.Sizes[0]):
                    self.grd[i] = self.Fijk ([i])
            elif self.dim == 2:
                for i in range(self.Sizes[0]):
                    for j in range(self.Sizes[1]):
                        self.grd[i,j] = self.Fijk ([i,j])


    def DERIV2_mine(self, d0, d1, argxy):   
        if self.SymbolDiffer:
            return self.Hessian[d0][d1](self.Node_to_Norm_or_Real(argxy)) 
        
        else:
            arg = np.array(self.node_to_normalized(argxy))  
            step = [a.step for a in self.A] 
            for i, el in enumerate(step):
                if i in self.grd_axis:
                    if self.grd_norm == "Node":
                        step[i] = 1
                    elif self.grd_norm == "Norm01":
                        step[i] = el / self.A[i].ma_mi
                if i in self.smb_axis:
                    if self.grd_norm == "Node":
                        step[i] = 1
                    elif self.grd_norm == "Norm01":
                        step[i] = el / self.A[i].ma_mi

            if d0==0 and d1==0 :
                if self.dim != 1:  step[1] = 0
                ret =  ( self.interpoleNode(arg-step) -2*self.interpoleNode(arg) +self.interpoleNode(arg+step) )/step[0]**2
            elif d0==1 and d1==1 :
                step[0] = 0
                ret =  ( self.interpoleNode(arg-step) -2*self.interpoleNode(arg) +self.interpoleNode(arg+step) )/step[1]**2
            else:      #if d0==0 and d1==1 :
                step1 = [step[0],-step[1]]
                ret =  ( self.interpoleNode(arg+step)+self.interpoleNode(arg-step)-self.interpoleNode(arg+step1)-self.interpoleNode(arg-step1) ) \
                    / step[0] / step[1] * .25
    
            return ret
            
    
    def F_mine ( self, ArS_real ) :  #   real args
        print("Вызывается F_my!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1")
        ret = self.interpoleNode( self.real_to_normalized( ArS_real) ) + self.V.avr    # gr = self.grd
        return ret

    # def Node_to_Norm_or_Real_mine (self, nodeArgs):
    #     #realArgs = self.Node_to_Real (nodeArgs)      
    #     #return self.Real_to_Norm_or_Real (realArgs)
    #     if self.ArgNorm :  return  [nodeArgs[i]*a.step / a.ma_mi   for i, a in enumerate(self.A)]
    #     else            :  return  [nodeArgs[i]*a.step + a.min     for i, a in enumerate(self.A)]
    
# #Из сеточных вариантов
#     def Ftbl ( self, n ) :
#         if self.type[0] == 'g':      # 2407
#             argsNode = [(a.dat[n]-a.min)/a.step   for a in self.A]
#  #       print (argsNode)
#             return self.interpolNode (argsNode)
        
#     def F ( self, ArS_real ) :
#           if self.dim == 0:                           #  31
#                 if SvF.Use_var: return self.var
#                 else          : return self.grd
#  #         print (self.ArgNorm)
#   #        1/0
#    #       if self.ArgNorm :  ArS_real = self.Norm01_to_Real(ArS_real)
#           if   len(ArS_real) > self.dim :
#               if ArS_real[self.dim] == 'N' :
#                   ArS_real = self.Norm01_to_Real(ArS_real)
#           elif SvF.F_Arg_Type == 'N' :  ArS_real = self.Norm01_to_Real(ArS_real)    #  заплатка для ArgNorm для fNi_fon(X,Y) символ функции   Ni(X,Y)  = Ni_fon(X,Y) + fon

#           ar =  self.Real_to_Node ( ArS_real )
#  #         print ('N', ar)
#           return self.interpolNode(ar)
# #          if ( self.type[0] == 'g' ) and self.dim==1 :       # 2407
#  #           return self.interpol ( 1, ar[0] )
#   #        elif ( self.type[0] == 'g' ) and self.dim==2 :       # 2407
#    #         return self.interpol ( 2, ar[0], ar[1] )
#     #      elif self.type[0] == 'g' and self.dim==3 :
#      #         return self.interpol ( 3, ar[0], ar[1], ar[2] )

# Вид функции для F и Ftbl должен быть: проверка на символьность - если да, то оставляю как сейчас. Если нет, то вызов interpolate. 
# Нужно лишь определиться с нормировкой и тем влияет ли она на что-то. Когда всё возврашается в обратный масштаб?

# Стоит отдельно нормировать сеточные и символьные. Это можно делать двумя способами: Nodes и [0,1]