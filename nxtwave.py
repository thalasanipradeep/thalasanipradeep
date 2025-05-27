import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class MobiusStrip:
    def _init_(self, R, W, resolution):
        if not all(isinstance(arg, (int, float)) and arg > 0 for arg in [R, W]):
            raise ValueError("R and W must be positive numbers.")
        if not isinstance(resolution, int) or resolution <= 0:
            raise ValueError("Resolution must be a positive integer.")
        self.R = R
        self.W = W
        self.resolution = resolution

        self.u = np.linspace(0, 2 * np.pi, self.resolution)
        self.v = np.linspace(-self.W / 2, self.W / 2, self.resolution)
        self.U, self.V = np.meshgrid(self.u, self.v)

        self._compute_3d_meshgrid()
        self._compute_surface_area()
        self._compute_edge_length()

    def _compute_3d_meshgrid(self):
        
        self.X = (self.R + self.V * np.cos(self.U / 2)) * np.cos(self.U)
        self.Y = (self.R + self.V * np.cos(self.U / 2)) * np.sin(self.U)
        self.Z = self.V * np.sin(self.U / 2)

    def _compute_surface_area(self):
        
        dXdu = - (self.R + self.V * np.cos(self.U / 2)) * np.sin(self.U) - \
               self.V / 2 * np.sin(self.U / 2) * np.cos(self.U)
        dYdu = (self.R + self.V * np.cos(self.U / 2)) * np.cos(self.U) - \
               self.V / 2 * np.sin(self.U / 2) * np.sin(self.U)
        dZdu = self.V / 2 * np.cos(self.U / 2)

        
        dXdv = np.cos(self.U / 2) * np.cos(self.U)
        dYdv = np.cos(self.U / 2) * np.sin(self.U)
        dZdv = np.sin(self.U / 2)

        E = dXdu*2 + dYdu + dZdu*2
        F = dXdu*dXdv + dYdu*dYdv + dZdu*dZdv
        G = dXdv*2 + dYdv + dZdv*2

        dA = np.sqrt(E * G - F**2)
        du_step = self.u[1] - self.u[0]
        dv_step = self.v[1] - self.v[0]
        self.surface_area = np.sum(dA) * du_step * dv_step

    def _compute_edge_length(self):
    
        V_edge_top = np.full_like(self.u, self.W / 2)
        X_edge_top = (self.R + V_edge_top * np.cos(self.u / 2)) * np.cos(self.u)
        Y_edge_top = (self.R + V_edge_top * np.cos(self.u / 2)) * np.sin(self.u)
        Z_edge_top = V_edge_top * np.sin(self.u / 2)
        V_edge_bottom = np.full_like(self.u, -self.W / 2)
        X_edge_bottom = (self.R + V_edge_bottom * np.cos(self.u / 2)) * np.cos(self.u)
        Y_edge_bottom = (self.R + V_edge_bottom * np.cos(self.u / 2)) * np.sin(self.u)
        Z_edge_bottom = V_edge_bottom * np.sin(self.u / 2)

        dX = np.diff(X_edge_top)
        dY = np.diff(Y_edge_top)
        dZ = np.diff(Z_edge_top)
        self.edge_length = np.sum(np.sqrt(dX*2 + dY + dZ*2))

    def get_mesh(self):
        
        return self.X, self.Y, self.Z

    def get_surface_area(self):
        
        return self.surface_area
    def get_edge_length(self):
        
        return self.edge_length

    def plot_mobius_strip(self):
        
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection='3d')

        ax.plot_surface(self.X, self.Y, self.Z, color='skyblue', edgecolor='none', alpha=0.8)
        V_edge_top = np.full_like(self.u, self.W / 2)
        X_edge_top = (self.R + V_edge_top * np.cos(self.u / 2)) * np.cos(self.u)
        Y_edge_top = (self.R + V_edge_top * np.cos(self.u / 2)) * np.sin(self.u)
        Z_edge_top = V_edge_top * np.sin(self.u / 2)

        V_edge_bottom = np.full_like(self.u, -self.W / 2)
        X_edge_bottom = (self.R + V_edge_bottom * np.cos(self.u / 2)) * np.cos(self.u)
        Y_edge_bottom = (self.R + V_edge_bottom * np.cos(self.u / 2)) * np.sin(self.u)
        Z_edge_bottom = V_edge_bottom * np.sin(self.u / 2)

        ax.plot(X_edge_top, Y_edge_top, Z_edge_top, color='red', linewidth=3, label='Edge (v=W/2)')
        ax.plot(X_edge_bottom, Y_edge_bottom, Z_edge_bottom, color='green', linewidth=3, linestyle='--', label='Edge (v=-W/2)')
        ax.set_title(f'Mobius Strip (R={self.R}, W={self.W})')
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.legend()
        plt.show()

    if __name__== "__main__":
      
      R_val = 2.0
      W_val = 1.0
      res_val = 100

    try:
        mobius = MobiusStrip(R=R_val, W=W_val, resolution=res_val)
        X, Y, Z = mobius.get_mesh()
        surface_area = mobius.get_surface_area()
        edge_length = mobius.get_edge_length()

        print(f"Mobius Strip Parameters: R={R_val}, W={W_val}, Resolution={res_val}")
        print(f"Approx. Surface Area: {surface_area:.4f}")
        print(f"Approx. Edge Length: {edge_length:.4f}")

        mobius.plot_mobius_strip()

        print("\n--- Another Mobius Strip Example ---")
        mobius_2 = MobiusStrip(R=3.0, W=0.8, resolution=80)
        print(f"Mobius Strip Parameters: R={mobius_2.R}, W={mobius_2.W}, Resolution={mobius_2.resolution}")
        print(f"Approx. Surface Area: {mobius_2.get_surface_area():.4f}")
        print(f"Approx. Edge Length: {mobius_2.get_edge_length():.4f}")
        mobius_2.plot_mobius_strip()

    except ValueError as e:
        print(f"Error: {e}")