from Lego import *
import math as m
from ModelFiles import to_logOut
import itertools

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

class MixedFun (BaseFun) :
    """
    Синтаксис: f ( X, V ) = [SPWL,Mesh,Polynome(6)];
    """

    def __init__(self, Vname='', As=[], discr_str=None, param=False, Finitialize=0, DataReadFrom='', Data=[], Type='', Domain=None, ReadFrom=''):
        Fun.__init__(self, Vname, As, param, -1, Finitialize, DataReadFrom, Data, Type, Domain, ReadFrom)

        self.discr_str = discr_str
        self.discr_list = dict()
        self.smb_discr_list = dict()
        self.grd_discr_list = dict()
        self.is_smb = True

        self.process_discr_string()
        self.print_info()
        # self.discr_types = discr_types
        self.dim = len(self.discr_list)

        self.discr_methods = { 
            'SPWL': self._SPWL_discr}
        
        # self.discr_methods = { 
        #     'SPWL': self._SPWL_discr,
        #     'Polynome': self._polynomial_discr,
        #     'Fourier': self._fourier_discr,
        #     'SOS2': self._SOS2_discr}
    
    def process_discr_string(self):
        # пусть на вход входит "[SPWL,Mesh,Polynome(6)]"
        discr_list = simple_split(self.discr_str[1:-1])
        for i, discr_method in enumerate(discr_list):
            if discr_method.startswith('SPWL') == 1 or self.discr_methods.startswith('Mesh') == 0:
                self.grd_discr_list[str(i)] = discr_method
            else:
                self.smb_discr_list[str(i)] = discr_method

        self.discr_list = discr_list
        if len(self.smb_discr_list) == 0: self.is_smb = False

        return 0
    
    def print_info(self):
        print("Function configuration:")
        print(f"    grid types = {self.grd_discr_list}")
        print(f"    symbolic types = {self.smb_discr_list}")

    
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
    

    def interpolNode ( self, argNode, lev=None):
        # изначально и вправду подаётся lev=None
        if lev is None :
            #print(self.A[0].Val, self.A[1].Val)
            #print(self.A[0].dat)
            lev = len(self.grd_discr_list)
            argNode1 = copy(argNode)
            argNode1 = [(argNode1[i] * self.A[i].step + self.A[i].min) for i in range(len(argNode))] # плохая строчка
            indices = []
        else:
            indices = copy(argNode)
            argNode1 = copy(argNode)

        if lev > 0 :
            print(f"lev={lev}", f"Node={argNode}", f"Node1={argNode1}")

            method_name = list(self.grd_discr_list.values())[lev-1] #[0]
            discr_func = self.discr_methods.get(method_name)
            axis_num = int(list(self.grd_discr_list.keys())[lev-1])

            if discr_func is None:
                raise ValueError(
                    f"Неизвестный метод дискретизации '{method_name}'. "
                    f"Доступные: {list(self.discr_methods.keys())}"
            )
            print(f"На этапе {lev} вызывается {discr_func.__name__}") # debug
            #print(f"текущий набор индексов {indices} и уровень {lev}")


            ret = discr_func(self.A[axis_num], [self.interpolNode([i] + indices, lev-1) 
                                           for i in self.A[axis_num].NodS], argNode1[axis_num])
            
            # print(f"self.A[lev-1] = {self.A[lev-1].Val}")
            # print(f"[self.interpolNode([i] + indices, lev-1) \
            #                                for i in self.A[lev-1].NodS] = {[self.interpolNode([i] + indices, lev-1) \
            #                                for i in self.A[lev-1].NodS]}")
            
            return ret
        
        else :
            print(f"lev={lev}", f"Node={argNode}", f"Node1={argNode1}")
            return self.fetch_value(argNode1)

    def fetch_value(self, argNode):
        if self.is_smb == False:
            return self.gr[tuple(argNode)]
        
        else:
            smb_axis = [int(_) for _ in self.smb_discr_list.keys()]
            needed_cords = [argNode[_] for _ in smb_axis]
            print(f"smb_axis = {smb_axis}")
            self.smbF(needed_cords)
        

#     def Ftbl ( self, n ) :
#     #     Args = [a.dat[n] for a in self.A]
#     # #    if self.ArgNorm:  SvF.F_Arg_Type = 'N'     #  заплатка для ArgNorm для fNi_fon(X,Y) символ функции   Ni(X,Y)  = Ni_fon(X,Y) + fon
#     #     ret = self.smbF( self.Real_to_Norm_or_Real( Args) )
#     ret = interp
#      #   SvF.F_Arg_Type = ''
#         return ret
# #        return self.smbF(self.Real_to_Norm_or_Real(Args))



    def F ( self, ArS_real ) :  # первая функция, вызываемая для поиска значения
        #ret = self.smbF( self.Real_to_Norm_or_Real( ArS_real) ) + self.V.avr    # gr = self.grd
        ret = self.interpolNode(ArS_real)

        return ret
    
    def Ftbl ( self, n ) :  # первая функция, вызываемая для поиска значения
        #ret = self.smbF( self.Real_to_Norm_or_Real( ArS_real) ) + self.V.avr    # gr = self.grd
        Args = [a.dat[n] for a in self.A]
        ret = self.interpolNode(Args)

        return ret