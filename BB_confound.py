import pandas as pd
import statsmodels.api as sm
from ConfoundingInterval import ConfoundingInterval
from time import time

class Node:
    def __init__(self, x_res, y_res, z_res, I_w, I_z, beta = None):
        self.beta = beta
        self.x_res = x_res
        self.y_res = y_res
        self.z_res = z_res
        self.I_w = I_w
        self.I_z = I_z

    def compute_CI(self, r2x, r2y):
        CI = ConfoundingInterval(self.y_res.corr(self.x_res), self.y_res.std(), self.x_res.std(), [0,0,-1], [r2x, r2y] + [1])
        CI.optimize('AOK')
        return [CI.fx_min, CI.fx_max]
    
class BB_confound:
    def __init__(self, x, y, s):
        self.start_time = time()
        self.x = x
        self.y = y
        self.s = s
        self.nodes_visited = 0
        self.temp_min = self.temp_max = sm.OLS(self.y, sm.add_constant(self.x)).fit().params.iloc[1]
        self.min_covariates = self.max_covariates = []
        self.queue = [Node(x - x.mean(), y - y.mean(), self.calc_z_res([], s.columns), [], s.columns)]

        self.bounds = self.get_bounds()

    def increment_nodes_visited(self):
        self.nodes_visited = self.nodes_visited + 1

    def calc_beta(self, node):
        if node.beta is None:
            node.beta = sm.OLS(node.y_res, sm.add_constant(node.x_res)).fit().params.iloc[1]
            if node.beta < self.temp_min:
                self.temp_min = node.beta
                self.min_covariates = node.I_w
            if node.beta > self.temp_max:
                self.temp_max = node.beta
                self.max_covariates = node.I_w
    
    def calc_r2_bounds(self,node):
        if len(node.I_z) > 0:
            return sm.OLS(node.x_res, sm.add_constant(node.z_res)).fit().rsquared, sm.OLS(node.y_res, sm.add_constant(node.z_res)).fit().rsquared
        return 0, 0
    
    def calc_zstar(self, node):
        return abs(node.z_res.corrwith(node.x_res, axis = 0) * node.z_res.corrwith(node.y_res, axis = 0)).idxmax()
    
    def calc_z_res(self, I_w, I_z):
        z_res = pd.DataFrame()
        if len(I_w) == 0:
            return self.s - self.s.mean()
        elif len(I_z) == 0:
            return None
        else:
            for z in I_z:
                z_res[z] = sm.OLS(self.s[z], sm.add_constant(self.s[I_w])).fit().resid
            return z_res
    
    def add_child_nodes(self, node):
        zstar = self.calc_zstar(node)
        self.queue.append(Node(node.x_res,
                               node.y_res,
                               node.z_res.drop(zstar, axis = 1),
                               node.I_w,
                               [z for z in node.I_z if z != zstar],
                               beta = node.beta))
        self.queue.append(Node(sm.OLS(self.x, sm.add_constant(self.s[node.I_w + [zstar]])).fit().resid,
                               sm.OLS(self.y, sm.add_constant(self.s[node.I_w + [zstar]])).fit().resid,
                               self.calc_z_res(node.I_w + [zstar], [z for z in node.I_z if z != zstar]),
                               node.I_w + [zstar],
                               [z for z in node.I_z if z != zstar]
        ))

    def get_bounds(self):
        while len(self.queue) > 0:
            self.increment_nodes_visited()
            current = self.queue.pop(0)
            self.calc_beta(current)
            r2x, r2y = self.calc_r2_bounds(current)
            lower, upper = current.compute_CI(r2x, r2y)
            if (lower < self.temp_min or upper > self.temp_max) and len(current.I_z) > 0:
                self.add_child_nodes(current)
        
        return [self.temp_min, self.temp_max]
    
    def __str__(self):
        output = f"Nodes visited: {self.nodes_visited} of {2**(len(self.s.columns) + 1) - 1} ({round(self.nodes_visited / (2**(len(self.s.columns) + 1) - 1) * 100, 2)}%)\n"
        output += f"Time elapsed:  {round(time() - self.start_time, 5)} seconds\n"
        output += f"Bounds: {round(self.temp_min, 6)}, {round(self.temp_max, 6)}\n\n"
        output += f"Min covariates: {self.min_covariates}\n"
        output += f"Max covariates: {self.max_covariates}"
        return output
