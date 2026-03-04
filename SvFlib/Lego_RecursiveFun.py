from Lego import *
import math as m
from ModelFiles import to_logOut
import itertools

class RecursiveFun (Fun) :
    """
    Класс, осуществляющий дискретизацию функции произвольного числа переменных
    рекурсивным проходом по осям. При этом для каждой оси задаётся свой тип 'одномерной'
    дискретизации из перечисленных: 
    1) SPWL - сглаженные кусочно-линейные функции
    2) pol - полиномы положительных степеней
    3) trig - частичные суммы ряда Фурье
    """

    def __init__(self, Vname='', As=[], discr_types=None, param=False, Finitialize=0, DataReadFrom='', Data=[], Type='', Domain=None, ReadFrom=''):
        Fun.__init__(self, Vname, As, param, -1, Finitialize, DataReadFrom, Data, Type, Domain, ReadFrom)

        self.discr_types = discr_types
        self.dim = len(discr_types)
        # with open('recursive_init.log', 'a') as f:
        #     f.write(f"\n{'='*60}\n")
        #     f.write(f"CONSTRUCTOR RecursiveFun called\n")
        #     f.write(f"  Vname: '{Vname}'\n")
        #     f.write(f"  As: {As}\n")
        #     f.write(f"  discr_types: {discr_types}\n")

        self.discr_methods = { 
            'SPWL': self._SPWL_discr,
            'Polynome': self._polynomial_discr,
            'Fourier': self._fourier_discr,
            'SOS2': self._SOS2_discr}


    @staticmethod
    def _CHKS_smoothing(z, e): # CHKS smoothing
        return 0.5 * (z + py.sqrt(z**2 + e))
    

    # def _SPWL_discr(self, x_points, y_points, x_eval, epsilon=0.1):
    #     step = x_points[1] - x_points[0] # для равномерной сетки

    #     derivatives = (y_points[1:] - y_points[:-1]) / step
    #     A_list = np.concatenate([[0], derivatives])

    #     result_start = y_points[0] + (A_list[1] - A_list[0]) * (x_eval - x_points[0]) 
        
    #     result_add = 0
    #     for i in range(2, len(y_points)):
    #         result_add += (A_list[i] - A_list[i-1]) * self._CHKS_smoothing(x_eval - x_points[i-1], epsilon)

    #     return result_start + result_add

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
        # print(x_points.Ub, len(y_points))

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
    def _SOS2_discr(x_points, y_points, x_eval):
        # Находим интервал
        print(x_points.Val, x_eval)
        i = x_points.ValToInd(x_eval)
        
        # Точное попадание в узел
        # if np.isclose(x_eval, x_points[i]):
        #     return y_points[i] #, [i], [1.0]
        
        if i <= 0: 
            return y_points[0]
    
        if i >= x_points.Ub: 
            return y_points[x_points.Ub]
        
        # Вычисляем веса для двух соседних точек (SOS2)
        x_i, x_i1 = x_points.Val[i], x_points.Val[i+1]
        d_total = x_points.step
        
        lambda_i = (x_i1 - x_eval) / d_total  # вес для i-й точки
        lambda_i1 = (x_eval - x_i) / d_total  # вес для (i+1)-й точки
        
        # Значение интерполяции
        value = lambda_i * y_points[i] + lambda_i1 * y_points[i+1]
        # to_logOut("Bruh")
        
        return value
    

    @staticmethod
    def _polynomial_discr(x_points: np.ndarray, y_points: np.ndarray, x_eval):
        n = x_points.Ub + 1
        
        denominators = [0 for _ in range(n)]
        for i in range(n):
            denominator = 1.0
            for j in range(n):
                if i != j:
                    denominator *= (x_points.Val[i] - x_points.Val[j])
            denominators[i] = denominator

        result = 0
        
        for i in range(n):
            # Вычисляем числитель для L_i(x)
            numerator = 1
            for j in range(n):
                if i != j:
                    numerator *= (x_eval - x_points.Val[j])
            
            # Базисный полином Лагранжа
            li = numerator / denominators[i]
            
            # Добавляем вклад узла
            result += y_points[i] * li
        
        return result

    @staticmethod
    def _fourier_discr(x_points: np.ndarray, y_points: np.ndarray, x_eval):
        """
        Интерполяция тригонометрический многочленом через ДПФ (дискретное преобразование Фурье). 
        Для хоть какой-то оптимизации используются векторные вычисления массивов numpy.
        """
        N = len(x_points)
        L = x_points[-1] - x_points[0]
        a = x_points[0]
        n = m.ceil((N + 1) / 2)
        
        # 1. Базовый угол для каждой точки
        base_angle = (2 * np.pi / L) * (x_points - a)  # размер N
        
        # 2. Создаём матрицу k для всех гармоник и всех точек
        # k_values: [0, 1, 2, ..., n-1] - вертикальный вектор
        k_values = np.arange(n).reshape(-1, 1)  # размер (n, 1)
        
        # 3. Вычисляем ВСЕ углы сразу: k * base_angle
        # Используем broadcasting: (n,1) * (1,N) → (n,N)
        all_angles = k_values @ base_angle.reshape(1, -1)  # матрица n×N
        
        # 4. Вычисляем ВСЕ косинусы и синусы сразу
        cos_matrix = np.cos(all_angles)  # матрица n×N
        sin_matrix = np.sin(all_angles)  # матрица n×N
        
        # 5. Умножаем на y_points и суммируем по точкам (ось 1)
        a_coeff = (2/N) * (cos_matrix @ y_points)  # вектор длины n
        b_coeff = (2/N) * (sin_matrix @ y_points)  # вектор длины n
        
        # 6. Поправки
        a_coeff[0] = a_coeff[0] / 2  # по формуле a_0, а не a_0/2
        b_coeff[0] = 0.0  # b_0 всегда 0
        
        if N % 2 == 0:  # чётное N
            b_coeff[-1] = 0.0  # b_n = 0

        new_coeff = (2 * np.pi / L) * (x_eval - a)    
        new_angles = (np.arange(n) * new_coeff)
        new_cos_vector = np.cos(new_angles)  # матрица n×N
        new_sin_vector = np.sin(new_angles)  # матрица n×N
        
        result = np.dot(new_sin_vector, b_coeff) + np.dot(new_cos_vector, a_coeff)
        
        # 9. Дополнительная поправка для чётного N
        if N % 2 == 0:
            k_last = n - 1
            # Вычитаем полный вклад и добавляем половинный
            correction = -0.5 * a_coeff[-1] * new_cos_vector[k_last]
            result += correction
        
        return result
    
    # NUMPY ВЕРСИЯ МЕТОДА
    # def interpolNode(self, argNode: list, lev=1):
    #     """
    #     Основной метод класса, считающий по переданной сетке значений значение интерполирующей функции в точке.
    #     Методом интерполяции является циклический вызов одномерных методов интерполяции, разных для
    #     каждой оси. 

    #     Parameters
    #     ----------
    #     point : list
    #         список с координатами точки, в которой рассчитывается значение

    #     Returns
    #     -------
    #     float
    #     Значение интерполирующей функции в указанной точке
    #     """        
        
    #     if len(argNode) != self.dim:
    #         raise ValueError(
    #             f"Число координат точки не соответсвует числу аргументов функции:  "
    #             f"{len(argNode)} вместо {self.dim}")
        
    #     if self.param or not SvF.Use_var:
    #         gr = self.grd
    #     else:
    #         gr = self.var

    #     for axes_num in range(self.dim):
    #         # Подготовка к проходу по оси
    #         print(f"Начинается проход по оси {axes_num}")
    #         if axes_num == 0:
    #             prev_data = gr

    #         method_name = self.discr_types[axes_num][0]
    #         discr_func = self.discr_methods.get(method_name)

    #         if discr_func is None:
    #             raise ValueError(
    #                 f"Неизвестный метод дискретизации '{method_name}'. "
    #                 f"Доступные: {list(self.discr_methods.keys())}"
    #         )
    #         print(f"На этапе {axes_num} вызывается {discr_func.__name__}") # debug

    #         other_sizes = prev_data.shape[1:] 
    #         new_data = np.zeros(other_sizes)

    #         # цикл интерполирования функции на следующий шаг
    #         for other_idx in np.ndindex(*other_sizes):
    #             slices = [slice(None)] 
    #             slices.extend([idx for idx in other_idx])
    #             x_slice = prev_data[tuple(slices)]
    #             print(f"индекс={other_idx}, слайс={x_slice}")

    #             if method_name == "SPWL" and len(self.discr_types[axes_num]) == 2:
    #                 print("epsilon")
    #                 custom_epsilon = self.discr_types[axes_num][1]
    #                 new_data[other_idx] = discr_func(self.A[axes_num], x_slice, argNode[axes_num], custom_epsilon)
    #             else:
    #                 new_data[other_idx] = discr_func(self.A[axes_num], x_slice, argNode[axes_num])

    #         prev_data = new_data
            

        return float(prev_data)

    # def interpolNode(self, argNode, lev=1):
    #     """
    #     Чисто Pyomo-версия интерполяции.
    #     Строит символическое выражение без numpy.
    #     """

    #     if len(argNode) != self.dim:
    #         raise ValueError(
    #             f"Число координат точки не соответствует размерности: "
    #             f"{len(argNode)} вместо {self.dim}"
    #         )

    #     # источник данных — только Pyomo
    #     gr = self.var

    #     # выбираем методы дискретизации по осям
    #     discr_funcs = []
    #     for axis in range(self.dim):
    #         method_name = self.discr_types[axis][0]
    #         discr_func = self.discr_methods.get(method_name)

    #         if discr_func is None:
    #             raise ValueError(
    #                 f"Неизвестный метод дискретизации '{method_name}'"
    #             )

    #         discr_funcs.append(discr_func)

    #     # --- основная идея ---
    #     # sum_{i1,i2,...,id} gr[i1,i2,...,id] *
    #     #   prod_k w_k(i_k, argNode[k])

    #     expr = 0

    #     # индексы всех узлов сетки
    #     index_sets = [self.A[k].NodSm for k in range(self.dim)]

    #     for idx in itertools.product(*index_sets):
    #         weight = 1

    #         for k in range(self.dim):
    #             discr_func = discr_funcs[k]

    #             if (
    #                 self.discr_types[k][0] == "SPWL"
    #                 and len(self.discr_types[k]) == 2
    #             ):
    #                 eps = self.discr_types[k][1]
    #                 weight *= discr_func(
    #                     self.A[k], idx[k], argNode[k], eps
    #                 )
    #             else:
    #                 weight *= discr_func(
    #                     self.A[k], idx[k], argNode[k]
    #                 )

    #         expr += gr[idx] * weight

    #     return expr

################################################################
    # def interpolNode(self, argNode, lev=None ):
    #     """
    #     Рекурсивно строит Pyomo выражение для интерполяции.
    #     """
    #     if self.param or not SvF.Use_var:
    #         gr = self.grd  # NumPy массив
    #     else:
    #         gr = self.var  # Pyomo IndexedVar

    #     def build_expr_for_axis(axis_idx, indices):
    #         """
    #         Рекурсивная функция построения выражения.
            
    #         Args:
    #             axis_idx: индекс текущей оси (0, 1, ...)
    #             indices: текущие индексы по уже обработанным осям
    #         """
    #         if axis_idx >= self.dim:
    #             # Все оси обработаны - возвращаем переменную Pyomo
    #             # Преобразуем индексы в tuple для доступа к Pyomo переменной
    #             index_tuple = tuple(int(idx) for idx in indices)
    #             return gr[index_tuple]
            
    #         method_name = self.discr_types[axis_idx][0]
    #         discr_func = self.discr_methods.get(method_name)
            
    #         if discr_func is None:
    #             raise ValueError(f"Неизвестный метод дискретизации '{method_name}'")
            
    #         # Определяем количество точек по текущей оси
    #         axis_points = self.A[axis_idx].Ub + 1
            
    #         if method_name == "SPWL" and len(self.discr_types[axis_idx]) == 2:
    #             custom_epsilon = self.discr_types[axis_idx][1]
    #             return discr_func(
    #                 self.A[axis_idx],
    #                 [build_expr_for_axis(axis_idx + 1, indices + [i]) 
    #                 for i in range(axis_points)],
    #                 argNode[axis_idx],
    #                 custom_epsilon
    #             )
    #         else:
    #             return discr_func(
    #                 self.A[axis_idx],
    #                 [build_expr_for_axis(axis_idx + 1, indices + [i]) 
    #                 for i in range(axis_points)],
    #                 argNode[axis_idx]
    #             )
        
    #     # Запускаем рекурсию с первой оси
    #     result_expr = build_expr_for_axis(0, [])
    #     return result_expr



    
    # def interpolNode(self, argNode, lev=None):
    #     """
    #     Рекурсивно строит Pyomo выражение для интерполяции.
    #     Работает с последней оси к первой.
    #     """
    #     if self.param or not SvF.Use_var:
    #         gr = self.grd  # NumPy массив
    #     else:
    #         gr = self.var  # Pyomo IndexedVar

    #     def build_expr_for_axis(axis_idx, fixed_indices):
    #         """
    #         Рекурсивная функция построения выражения.
            
    #         Args:
    #             axis_idx: индекс текущей оси (начинаем с 0, но логика обратная)
    #             fixed_indices: уже фиксированные индексы по ПРЕДЫДУЩИМ осям
    #                         (те, что ближе к началу)
    #         """
    #         # Если это последняя ось (axis_idx == self.dim-1)
    #         # Мы должны собрать значения по этой оси для всех её точек
    #         if axis_idx == self.dim - 1:
    #             # Собираем Pyomo переменные для всех точек по последней оси
    #             values = []
    #             for i in range(self.A[axis_idx].Ub + 1):
    #                 # Собираем полный индекс: fixed_indices + [i]
    #                 full_index = tuple(list(fixed_indices) + [i])
    #                 values.append(gr[full_index])
                
    #             # Применяем интерполяцию по последней оси
    #             method_name = self.discr_types[axis_idx][0]
    #             discr_func = self.discr_methods.get(method_name)
                
    #             if discr_func is None:
    #                 raise ValueError(f"Неизвестный метод '{method_name}'")
                
    #             if method_name == "SPWL" and len(self.discr_types[axis_idx]) == 2:
    #                 custom_epsilon = self.discr_types[axis_idx][1]
    #                 return discr_func(
    #                     self.A[axis_idx],
    #                     values,
    #                     argNode[axis_idx],
    #                     custom_epsilon
    #                 )
    #             else:
    #                 return discr_func(
    #                     self.A[axis_idx],
    #                     values,
    #                     argNode[axis_idx]
    #                 )
            
    #         # Если не последняя ось, рекурсивно строим выражения
    #         # для каждой точки по текущей оси
    #         else:
    #             expressions = []
    #             for i in range(self.A[axis_idx].Ub + 1):
    #                 # Фиксируем текущую ось и идём глубже
    #                 expr = build_expr_for_axis(
    #                     axis_idx + 1, 
    #                     fixed_indices + [i]  # Добавляем индекс текущей оси
    #                 )
    #                 expressions.append(expr)
                
    #             # Теперь интерполируем по текущей оси
    #             method_name = self.discr_types[axis_idx][0]
    #             discr_func = self.discr_methods.get(method_name)
                
    #             if discr_func is None:
    #                 raise ValueError(f"Неизвестный метод '{method_name}'")
                
    #             if method_name == "SPWL" and len(self.discr_types[axis_idx]) == 2:
    #                 custom_epsilon = self.discr_types[axis_idx][1]
    #                 return discr_func(
    #                     self.A[axis_idx],
    #                     expressions,  # Список УЖЕ интерполированных выражений!
    #                     argNode[axis_idx],
    #                     custom_epsilon
    #                 )
    #             else:
    #                 return discr_func(
    #                     self.A[axis_idx],
    #                     expressions,
    #                     argNode[axis_idx]
    #                 )
        
    #     # Запускаем рекурсию с ПЕРВОЙ оси (индекс 0)
    #     # Но логика будет идти от последней к первой
    #     return build_expr_for_axis(0, [])
    

    def interpolNode ( self, argNode, lev=None):
        # изначально и вправду подаётся lev=None
        if self.param or not SvF.Use_var:
            gr = self.grd  # NumPy массив
        else:
            gr = self.var  # Pyomo IndexedVar

        if lev is None :
            #print(self.A[0].Val, self.A[1].Val)
            #print(self.A[0].dat)
            lev = len(argNode)
            argNode1 = copy(argNode)
            argNode1 = [(argNode1[i] * self.A[i].step + self.A[i].min) for i in range(len(argNode))]
            indices = []
        else:
            indices = copy(argNode)
            argNode1 = copy(argNode)
        if lev > 0 :
            print(f"lev={lev}", f"Node={argNode}", f"Node1={argNode1}")

            method_name = self.discr_types[lev-1] #[0]
            discr_func = self.discr_methods.get(method_name)

            if discr_func is None:
                raise ValueError(
                    f"Неизвестный метод дискретизации '{method_name}'. "
                    f"Доступные: {list(self.discr_methods.keys())}"
            )
            print(f"На этапе {lev} вызывается {discr_func.__name__}") # debug
            #print(f"текущий набор индексов {indices} и уровень {lev}")


            ret = discr_func(self.A[lev-1], [self.interpolNode([i] + indices, lev-1) 
                                           for i in self.A[lev-1].NodS], argNode1[lev-1])
            
            # print(f"self.A[lev-1] = {self.A[lev-1].Val}")
            # print(f"[self.interpolNode([i] + indices, lev-1) \
            #                                for i in self.A[lev-1].NodS] = {[self.interpolNode([i] + indices, lev-1) \
            #                                for i in self.A[lev-1].NodS]}")
            
            return ret
        else :
            #print(f"текущий набор индексов {argNode} и уровень {lev}")
            # print(argNode)
            print(f"lev={lev}", f"Node={argNode}", f"Node1={argNode1}")
            return self.gr[tuple(argNode1)]
        
    # def interpolNode(self, argNode, lev=None):
    #     if lev is None:
    #         lev = len(argNode)
        
    #     if lev > 0:
    #         # ВСЕГДА копируем argNode
    #         current_argNode = argNode.copy()
            
    #         method_name = self.discr_types[lev-1][0]
    #         discr_func = self.discr_methods.get(method_name)
            
    #         if discr_func is None:
    #             raise ValueError(f"Неизвестный метод '{method_name}'")
            
    #         # Собираем значения для текущей оси
    #         values = []
    #         for i in self.A[lev-1].NodS:
    #             # Фиксируем текущую ось
    #             new_argNode = current_argNode.copy()
    #             new_argNode[lev-1] = i
                
    #             # Рекурсивный вызов
    #             values.append(self.interpolNode(new_argNode, lev-1))
            
    #         # Применяем интерполяцию
    #         if method_name == "SPWL" and len(self.discr_types[lev-1]) == 2:
    #             return discr_func(
    #                 self.A[lev-1],
    #                 values,
    #                 current_argNode[lev-1],
    #                 self.discr_types[lev-1][1]
    #             )
    #         else:
    #             return discr_func(
    #                 self.A[lev-1],
    #                 values,
    #                 current_argNode[lev-1]
    #             )
    #     else:
    #         # Базовый случай
    #         return self.gr[tuple(argNode)]