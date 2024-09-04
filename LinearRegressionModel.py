import tensorflow as tf

class LinearRegressionModel(tf.keras.Model):
    def __init__(self, A, b, calc_resid = False, calc_rsquared = False):
        super(LinearRegressionModel, self).__init__()
        self.A = tf.pad(tf.cast(A, tf.float64), [[0,0],[1,0]], constant_values = 1)
        self.b = tf.cast(b, tf.float64)
        self.params = self.get_params()
        if calc_resid:
            self.resid = self.get_resid()
        else:
            self.resid = None
        if calc_resid:
            self.rsquared = self.get_rsquared()
        else:
            self.rsquared = None
        self.pred = None

    def call(self):
        if self.pred is None:
            self.pred = tf.matmul(self.A, self.params)
        return self.pred
    
    def get_rsquared(self):
        if self.pred is None:
            self.call()
        self.rsquared = 1 - tf.reduce_sum(tf.square(self.b - self.pred)) / tf.reduce_sum(tf.square(self.b - tf.reduce_mean(self.b)))
        return tf.squeeze(self.rsquared).numpy()
    
    def get_resid(self):
        if self.pred is None:
            self.call()
        self.resid = self.b - self.pred
        return self.resid
    
    def get_params(self, tol = 1e-8):
        n = self.A.get_shape()[1]
        x = tf.cast(tf.zeros([n,1]), self.A.dtype)
        r = self.b
        s = p = tf.matmul(tf.transpose(self.A), r)
        gamma = self.sq_norm2(s)
        for _ in range(5 * n):
            q = tf.matmul(self.A, p)
            alpha = gamma / self.sq_norm2(q)
            x = x + alpha * p
            r_new = r - alpha * q
            if self.sq_norm2(r) - self.sq_norm2(r_new) < tol:
                return x
            r = r_new
            s = tf.matmul(tf.transpose(self.A), r)
            gamma_new = self.sq_norm2(s)
            beta = gamma_new / gamma
            gamma = gamma_new
            p = s + beta * p
        self.params = x
        return x
    
    def sq_norm2(self, v):
        return tf.matmul(tf.transpose(v),v)
