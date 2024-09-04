from ConfoundingInterval import ConfoundingInterval
from time import time
import tensorflow as tf
import tensorflow_probability as tfp
from LinearRegressionModel import LinearRegressionModel

class Node:
    def __init__(self, x_res, y_res, z_res, I_w, I_z, beta = None):
        self.beta = beta
        self.x_res = x_res
        self.y_res = y_res
        self.z_res = z_res
        self.I_w = I_w
        self.I_z = I_z

    def compute_CI(self, r2x, r2y):
        CI = ConfoundingInterval(tf.squeeze(tfp.stats.correlation(self.x_res,self.y_res)).numpy(), tf.math.reduce_std(self.y_res).numpy(), tf.math.reduce_std(self.x_res).numpy(), [0,0,-1], [r2x,r2y,1])
        CI.optimize('AOK')
        return [CI.fx_min, CI.fx_max]
    
class BB_confound_tf:
    def __init__(self, x, y, s):
        self.start_time = time()
        self.col_names = s.columns
        self.x = tf.cast(tf.reshape(tf.convert_to_tensor(x), [len(x),1]), tf.float64)
        self.y = tf.cast(tf.reshape(tf.convert_to_tensor(y), [len(x),1]), tf.float64)
        self.s = tf.cast(tf.convert_to_tensor(s), tf.float64)
        self.nodes_visited = 0
        self.temp_min = self.temp_max = LinearRegressionModel(self.x,self.y).params[1,0].numpy()
        self.min_covariates = self.max_covariates = []
        self.queue = [Node(self.x - tf.reduce_mean(self.x),
                           self.y - tf.reduce_mean(self.y),
                           self.calc_z_res([], list(range(self.s.shape[1]))),
                           [],
                           list(range(self.s.shape[1]))
        )]

        self.bounds = self.get_bounds()

    def increment_nodes_visited(self):
        self.nodes_visited = self.nodes_visited + 1

    def calc_beta(self, node):
        if node.beta is None:
            node.beta = LinearRegressionModel(node.x_res,node.y_res).params[1,0].numpy()
            if node.beta < self.temp_min:
                self.temp_min = node.beta
                self.min_covariates = node.I_w
            if node.beta > self.temp_max:
                self.temp_max = node.beta
                self.max_covariates = node.I_w
    
    def calc_r2_bounds(self,node):
        if len(node.I_z) > 0:
            return LinearRegressionModel(node.z_res,node.x_res).get_rsquared(), LinearRegressionModel(node.z_res,node.y_res).get_rsquared()
        return 0, 0
    
    def calc_zstar(self, node):
        return tf.math.argmax(tf.math.abs(tf.map_fn(lambda col1: tfp.stats.correlation(col1,node.x_res[:,0],event_axis = None), tf.transpose(node.z_res)) *
                                          tf.map_fn(lambda col2: tfp.stats.correlation(col2,node.y_res[:,0],event_axis = None), tf.transpose(node.z_res)))).numpy()
    
    def calc_z_res(self, I_w, I_z):
        if len(I_w) == 0:
            return self.s - tf.reduce_mean(self.s)
        elif len(I_z) == 0:
            return None
        else:
            z_res = tf.transpose(tf.squeeze(tf.map_fn(lambda col: LinearRegressionModel(tf.gather(self.s, I_w, axis = 1), tf.reshape(col, [self.s.shape[0],1])).get_resid(), tf.transpose(tf.gather(self.s, I_z, axis = 1)))))
            if z_res.ndim == 1:
                z_res = tf.reshape(z_res, [z_res.shape[0],1])
            return z_res
    
    def add_child_nodes(self, node):
        zstar = self.calc_zstar(node)
        z_res_index = [z for z in range(len(node.I_z)) if z != zstar]
        I_z = [node.I_z[z] for z in z_res_index]
        # node1
        self.queue.append(Node(node.x_res,
                               node.y_res,
                               tf.gather(node.z_res, z_res_index, axis = 1) if len(I_z) > 0 else None,
                               node.I_w,
                               I_z,
                               beta = node.beta))
        # node2
        I_w = node.I_w + [node.I_z[zstar]]
        self.queue.append(Node(LinearRegressionModel(tf.gather(self.s, I_w, axis = 1), self.x).get_resid(),
                               LinearRegressionModel(tf.gather(self.s, I_w, axis = 1), self.y).get_resid(),
                               self.calc_z_res(I_w, I_z),
                               I_w,
                               I_z
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
        output = f"Nodes visited: {self.nodes_visited} of {2**(self.s.shape[1] + 1) - 1} ({round(self.nodes_visited / (2**(self.s.shape[1] + 1) - 1) * 100, 2)}%)\n"
        output += f"Time elapsed:  {round(time() - self.start_time, 5)} seconds\n"
        output += f"Bounds: {round(self.temp_min, 6)}, {round(self.temp_max, 6)}\n\n"
        output += f"Min covariates: {list(self.col_names[self.min_covariates])}\n"
        output += f"Max covariates: {list(self.col_names[self.max_covariates])}"
        return output
