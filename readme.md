# Insulating Layer PISO Solver

## 项目简介

该项目实现了一个使用 PISO 方法求解绝缘层内流体动力学和热传递的有限体积法 (FVM) 求解器。项目包括通过 Python 编写的 FVM 代码来模拟流体的速度场、压力场以及温度分布。所有的计算均基于径向网格，并生成相应的图像来展示计算结果。

## 文件结构

- [`MeshFVM.py`](command:_github.copilot.openRelativePath?%5B%7B%22scheme%22%3A%22file%22%2C%22authority%22%3A%22%22%2C%22path%22%3A%22%2Fd%3A%2FDesktop%2F%E4%B8%AD%E5%9B%BD%E5%BB%BA%E6%9D%90%2FCFD%E7%83%AD%E8%9C%A1%2FInsulating_layer_PISO%2FMeshFVM.py%22%2C%22query%22%3A%22%22%2C%22fragment%22%3A%22%22%7D%2C%221e544a3e-16c3-4e92-bcb9-b8326c54e279%22%5D "d:\Desktop\中国建材\CFD热蜡\Insulating_layer_PISO\MeshFVM.py")：用于生成径向网格并定义网格参数。
  - [`MeshStructure`](command:_github.copilot.openSymbolFromReferences?%5B%22%22%2C%5B%7B%22uri%22%3A%7B%22scheme%22%3A%22file%22%2C%22authority%22%3A%22%22%2C%22path%22%3A%22%2Fd%3A%2FDesktop%2F%E4%B8%AD%E5%9B%BD%E5%BB%BA%E6%9D%90%2FCFD%E7%83%AD%E8%9C%A1%2FInsulating_layer_PISO%2FMeshFVM.py%22%2C%22query%22%3A%22%22%2C%22fragment%22%3A%22%22%7D%2C%22pos%22%3A%7B%22line%22%3A2%2C%22character%22%3A6%7D%7D%5D%2C%221e544a3e-16c3-4e92-bcb9-b8326c54e279%22%5D "Go to definition") 类：定义网格结构。
  - [`MeshGenerator`](command:_github.copilot.openSymbolFromReferences?%5B%22%22%2C%5B%7B%22uri%22%3A%7B%22scheme%22%3A%22file%22%2C%22authority%22%3A%22%22%2C%22path%22%3A%22%2Fd%3A%2FDesktop%2F%E4%B8%AD%E5%9B%BD%E5%BB%BA%E6%9D%90%2FCFD%E7%83%AD%E8%9C%A1%2FInsulating_layer_PISO%2FMeshFVM.py%22%2C%22query%22%3A%22%22%2C%22fragment%22%3A%22%22%7D%2C%22pos%22%3A%7B%22line%22%3A12%2C%22character%22%3A6%7D%7D%2C%7B%22uri%22%3A%7B%22scheme%22%3A%22file%22%2C%22authority%22%3A%22%22%2C%22path%22%3A%22%2FD%3A%2FDesktop%2F%E4%B8%AD%E5%9B%BD%E5%BB%BA%E6%9D%90%2FCFD%E7%83%AD%E8%9C%A1%2FInsulating_layer_PISO%2FNS_FVM.py%22%2C%22query%22%3A%22%22%2C%22fragment%22%3A%22%22%7D%2C%22pos%22%3A%7B%22line%22%3A72%2C%22character%22%3A8%7D%7D%5D%2C%221e544a3e-16c3-4e92-bcb9-b8326c54e279%22%5D "Go to definition") 类：生成不同类型的网格，包括均匀网格、径向网格和非均匀网格。
- [`NS_FVM.py`](command:_github.copilot.openRelativePath?%5B%7B%22scheme%22%3A%22file%22%2C%22authority%22%3A%22%22%2C%22path%22%3A%22%2Fd%3A%2FDesktop%2F%E4%B8%AD%E5%9B%BD%E5%BB%BA%E6%9D%90%2FCFD%E7%83%AD%E8%9C%A1%2FInsulating_layer_PISO%2FNS_FVM.py%22%2C%22query%22%3A%22%22%2C%22fragment%22%3A%22%22%7D%2C%221e544a3e-16c3-4e92-bcb9-b8326c54e279%22%5D "d:\Desktop\中国建材\CFD热蜡\Insulating_layer_PISO\NS_FVM.py")：求解 Navier-Stokes 方程的核心代码，实现了速度、压力的计算。
  - 包含物理模型参数和数值模拟参数的设置。
  - 创建径向网格并进行坐标转换。
  - 初始化变量和边界条件。
  - 迭代求解 Navier-Stokes 方程并生成结果图像。
- [`solver.py`](command:_github.copilot.openRelativePath?%5B%7B%22scheme%22%3A%22file%22%2C%22authority%22%3A%22%22%2C%22path%22%3A%22%2Fd%3A%2FDesktop%2F%E4%B8%AD%E5%9B%BD%E5%BB%BA%E6%9D%90%2FCFD%E7%83%AD%E8%9C%A1%2FInsulating_layer_PISO%2Fsolver.py%22%2C%22query%22%3A%22%22%2C%22fragment%22%3A%22%22%7D%2C%221e544a3e-16c3-4e92-bcb9-b8326c54e279%22%5D "d:\Desktop\中国建材\CFD热蜡\Insulating_layer_PISO\solver.py")：调用 Navier-Stokes 求解器并处理温度场的求解。
  - 包含多个求解器函数，如 [`U_Solver`](command:_github.copilot.openSymbolFromReferences?%5B%22%22%2C%5B%7B%22uri%22%3A%7B%22scheme%22%3A%22file%22%2C%22authority%22%3A%22%22%2C%22path%22%3A%22%2Fd%3A%2FDesktop%2F%E4%B8%AD%E5%9B%BD%E5%BB%BA%E6%9D%90%2FCFD%E7%83%AD%E8%9C%A1%2FInsulating_layer_PISO%2Fsolver.py%22%2C%22query%22%3A%22%22%2C%22fragment%22%3A%22%22%7D%2C%22pos%22%3A%7B%22line%22%3A8%2C%22character%22%3A4%7D%7D%2C%7B%22uri%22%3A%7B%22scheme%22%3A%22file%22%2C%22authority%22%3A%22%22%2C%22path%22%3A%22%2FD%3A%2FDesktop%2F%E4%B8%AD%E5%9B%BD%E5%BB%BA%E6%9D%90%2FCFD%E7%83%AD%E8%9C%A1%2FInsulating_layer_PISO%2FNS_FVM.py%22%2C%22query%22%3A%22%22%2C%22fragment%22%3A%22%22%7D%2C%22pos%22%3A%7B%22line%22%3A194%2C%22character%22%3A41%7D%7D%5D%2C%221e544a3e-16c3-4e92-bcb9-b8326c54e279%22%5D "Go to definition")、[`V_Solver`](command:_github.copilot.openSymbolFromReferences?%5B%22%22%2C%5B%7B%22uri%22%3A%7B%22scheme%22%3A%22file%22%2C%22authority%22%3A%22%22%2C%22path%22%3A%22%2Fd%3A%2FDesktop%2F%E4%B8%AD%E5%9B%BD%E5%BB%BA%E6%9D%90%2FCFD%E7%83%AD%E8%9C%A1%2FInsulating_layer_PISO%2Fsolver.py%22%2C%22query%22%3A%22%22%2C%22fragment%22%3A%22%22%7D%2C%22pos%22%3A%7B%22line%22%3A110%2C%22character%22%3A4%7D%7D%2C%7B%22uri%22%3A%7B%22scheme%22%3A%22file%22%2C%22authority%22%3A%22%22%2C%22path%22%3A%22%2FD%3A%2FDesktop%2F%E4%B8%AD%E5%9B%BD%E5%BB%BA%E6%9D%90%2FCFD%E7%83%AD%E8%9C%A1%2FInsulating_layer_PISO%2FNS_FVM.py%22%2C%22query%22%3A%22%22%2C%22fragment%22%3A%22%22%7D%2C%22pos%22%3A%7B%22line%22%3A196%2C%22character%22%3A41%7D%7D%5D%2C%221e544a3e-16c3-4e92-bcb9-b8326c54e279%22%5D "Go to definition")、[`T_Solver`](command:_github.copilot.openSymbolFromReferences?%5B%22%22%2C%5B%7B%22uri%22%3A%7B%22scheme%22%3A%22file%22%2C%22authority%22%3A%22%22%2C%22path%22%3A%22%2Fd%3A%2FDesktop%2F%E4%B8%AD%E5%9B%BD%E5%BB%BA%E6%9D%90%2FCFD%E7%83%AD%E8%9C%A1%2FInsulating_layer_PISO%2Fsolver.py%22%2C%22query%22%3A%22%22%2C%22fragment%22%3A%22%22%7D%2C%22pos%22%3A%7B%22line%22%3A219%2C%22character%22%3A4%7D%7D%2C%7B%22uri%22%3A%7B%22scheme%22%3A%22file%22%2C%22authority%22%3A%22%22%2C%22path%22%3A%22%2FD%3A%2FDesktop%2F%E4%B8%AD%E5%9B%BD%E5%BB%BA%E6%9D%90%2FCFD%E7%83%AD%E8%9C%A1%2FInsulating_layer_PISO%2FNS_FVM.py%22%2C%22query%22%3A%22%22%2C%22fragment%22%3A%22%22%7D%2C%22pos%22%3A%7B%22line%22%3A214%2C%22character%22%3A34%7D%7D%5D%2C%221e544a3e-16c3-4e92-bcb9-b8326c54e279%22%5D "Go to definition") 和 [`P_Solver`](command:_github.copilot.openSymbolFromReferences?%5B%22%22%2C%5B%7B%22uri%22%3A%7B%22scheme%22%3A%22file%22%2C%22authority%22%3A%22%22%2C%22path%22%3A%22%2Fd%3A%2FDesktop%2F%E4%B8%AD%E5%9B%BD%E5%BB%BA%E6%9D%90%2FCFD%E7%83%AD%E8%9C%A1%2FInsulating_layer_PISO%2Fsolver.py%22%2C%22query%22%3A%22%22%2C%22fragment%22%3A%22%22%7D%2C%22pos%22%3A%7B%22line%22%3A358%2C%22character%22%3A4%7D%7D%2C%7B%22uri%22%3A%7B%22scheme%22%3A%22file%22%2C%22authority%22%3A%22%22%2C%22path%22%3A%22%2FD%3A%2FDesktop%2F%E4%B8%AD%E5%9B%BD%E5%BB%BA%E6%9D%90%2FCFD%E7%83%AD%E8%9C%A1%2FInsulating_layer_PISO%2FNS_FVM.py%22%2C%22query%22%3A%22%22%2C%22fragment%22%3A%22%22%7D%2C%22pos%22%3A%7B%22line%22%3A199%2C%22character%22%3A19%7D%7D%5D%2C%221e544a3e-16c3-4e92-bcb9-b8326c54e279%22%5D "Go to definition")。
  - 使用并行化和优化技术加速计算。
- [`Pressure.png`](command:_github.copilot.openRelativePath?%5B%7B%22scheme%22%3A%22file%22%2C%22authority%22%3A%22%22%2C%22path%22%3A%22%2Fd%3A%2FDesktop%2F%E4%B8%AD%E5%9B%BD%E5%BB%BA%E6%9D%90%2FCFD%E7%83%AD%E8%9C%A1%2FInsulating_layer_PISO%2FPressure.png%22%2C%22query%22%3A%22%22%2C%22fragment%22%3A%22%22%7D%2C%221e544a3e-16c3-4e92-bcb9-b8326c54e279%22%5D "d:\Desktop\中国建材\CFD热蜡\Insulating_layer_PISO\Pressure.png")：压力场分布图。
- [`RadialMesh.png`](command:_github.copilot.openRelativePath?%5B%7B%22scheme%22%3A%22file%22%2C%22authority%22%3A%22%22%2C%22path%22%3A%22%2Fd%3A%2FDesktop%2F%E4%B8%AD%E5%9B%BD%E5%BB%BA%E6%9D%90%2FCFD%E7%83%AD%E8%9C%A1%2FInsulating_layer_PISO%2FRadialMesh.png%22%2C%22query%22%3A%22%22%2C%22fragment%22%3A%22%22%7D%2C%221e544a3e-16c3-4e92-bcb9-b8326c54e279%22%5D "d:\Desktop\中国建材\CFD热蜡\Insulating_layer_PISO\RadialMesh.png")：径向网格的可视化。
- [`Temperature.png`](command:_github.copilot.openRelativePath?%5B%7B%22scheme%22%3A%22file%22%2C%22authority%22%3A%22%22%2C%22path%22%3A%22%2Fd%3A%2FDesktop%2F%E4%B8%AD%E5%9B%BD%E5%BB%BA%E6%9D%90%2FCFD%E7%83%AD%E8%9C%A1%2FInsulating_layer_PISO%2FTemperature.png%22%2C%22query%22%3A%22%22%2C%22fragment%22%3A%22%22%7D%2C%221e544a3e-16c3-4e92-bcb9-b8326c54e279%22%5D "d:\Desktop\中国建材\CFD热蜡\Insulating_layer_PISO\Temperature.png")：温度分布图。
- [`VelocityQuiver.png`](command:_github.copilot.openRelativePath?%5B%7B%22scheme%22%3A%22file%22%2C%22authority%22%3A%22%22%2C%22path%22%3A%22%2Fd%3A%2FDesktop%2F%E4%B8%AD%E5%9B%BD%E5%BB%BA%E6%9D%90%2FCFD%E7%83%AD%E8%9C%A1%2FInsulating_layer_PISO%2FVelocityQuiver.png%22%2C%22query%22%3A%22%22%2C%22fragment%22%3A%22%22%7D%2C%221e544a3e-16c3-4e92-bcb9-b8326c54e279%22%5D "d:\Desktop\中国建材\CFD热蜡\Insulating_layer_PISO\VelocityQuiver.png")：速度场矢量图。
- [`XVelocity.png`](command:_github.copilot.openRelativePath?%5B%7B%22scheme%22%3A%22file%22%2C%22authority%22%3A%22%22%2C%22path%22%3A%22%2Fd%3A%2FDesktop%2F%E4%B8%AD%E5%9B%BD%E5%BB%BA%E6%9D%90%2FCFD%E7%83%AD%E8%9C%A1%2FInsulating_layer_PISO%2FXVelocity.png%22%2C%22query%22%3A%22%22%2C%22fragment%22%3A%22%22%7D%2C%221e544a3e-16c3-4e92-bcb9-b8326c54e279%22%5D "d:\Desktop\中国建材\CFD热蜡\Insulating_layer_PISO\XVelocity.png")：X方向速度场分布图。
- [`YVelocity.png`](command:_github.copilot.openRelativePath?%5B%7B%22scheme%22%3A%22file%22%2C%22authority%22%3A%22%22%2C%22path%22%3A%22%2Fd%3A%2FDesktop%2F%E4%B8%AD%E5%9B%BD%E5%BB%BA%E6%9D%90%2FCFD%E7%83%AD%E8%9C%A1%2FInsulating_layer_PISO%2FYVelocity.png%22%2C%22query%22%3A%22%22%2C%22fragment%22%3A%22%22%7D%2C%221e544a3e-16c3-4e92-bcb9-b8326c54e279%22%5D "d:\Desktop\中国建材\CFD热蜡\Insulating_layer_PISO\YVelocity.png")：Y方向速度场分布图。
- [`MassResidual.png`](command:_github.copilot.openRelativePath?%5B%7B%22scheme%22%3A%22file%22%2C%22authority%22%3A%22%22%2C%22path%22%3A%22%2Fd%3A%2FDesktop%2F%E4%B8%AD%E5%9B%BD%E5%BB%BA%E6%9D%90%2FCFD%E7%83%AD%E8%9C%A1%2FInsulating_layer_PISO%2FMassResidual.png%22%2C%22query%22%3A%22%22%2C%22fragment%22%3A%22%22%7D%2C%221e544a3e-16c3-4e92-bcb9-b8326c54e279%22%5D "d:\Desktop\中国建材\CFD热蜡\Insulating_layer_PISO\MassResidual.png"): 质量残差图。
- [`XvelocityResidual.png`](command:_github.copilot.openRelativePath?%5B%7B%22scheme%22%3A%22file%22%2C%22authority%22%3A%22%22%2C%22path%22%3A%22%2Fd%3A%2FDesktop%2F%E4%B8%AD%E5%9B%BD%E5%BB%BA%E6%9D%90%2FCFD%E7%83%AD%E8%9C%A1%2FInsulating_layer_PISO%2FXvelocityResidual.png%22%2C%22query%22%3A%22%22%2C%22fragment%22%3A%22%22%7D%2C%221e544a3e-16c3-4e92-bcb9-b8326c54e279%22%5D "d:\Desktop\中国建材\CFD热蜡\Insulating_layer_PISO\XvelocityResidual.png"): X方向速度残差图。
- [`YvelocityResidual.png`](command:_github.copilot.openRelativePath?%5B%7B%22scheme%22%3A%22file%22%2C%22authority%22%3A%22%22%2C%22path%22%3A%22%2Fd%3A%2FDesktop%2F%E4%B8%AD%E5%9B%BD%E5%BB%BA%E6%9D%90%2FCFD%E7%83%AD%E8%9C%A1%2FInsulating_layer_PISO%2FYvelocityResidual.png%22%2C%22query%22%3A%22%22%2C%22fragment%22%3A%22%22%7D%2C%221e544a3e-16c3-4e92-bcb9-b8326c54e279%22%5D "d:\Desktop\中国建材\CFD热蜡\Insulating_layer_PISO\YvelocityResidual.png"): Y方向速度残差图。

## 依赖项

运行该项目需要以下依赖库：

```bash
pip install numpy matplotlib scipy numba
```

## 安装指南

1. 克隆该项目到本地：

```bash
git clone https://github.com/yourusername/insulating-layer-piso-solver.git
cd insulating-layer-piso-solver
```

2. 安装依赖库：

```bash
pip install numpy matplotlib scipy numba
```

## 使用说明

1. 运行 [`NS_FVM.py`](command:_github.copilot.openRelativePath?%5B%7B%22scheme%22%3A%22file%22%2C%22authority%22%3A%22%22%2C%22path%22%3A%22%2Fd%3A%2FDesktop%2F%E4%B8%AD%E5%9B%BD%E5%BB%BA%E6%9D%90%2FCFD%E7%83%AD%E8%9C%A1%2FInsulating_layer_PISO%2FNS_FVM.py%22%2C%22query%22%3A%22%22%2C%22fragment%22%3A%22%22%7D%2C%221e544a3e-16c3-4e92-bcb9-b8326c54e279%22%5D "d:\Desktop\中国建材\CFD热蜡\Insulating_layer_PISO\NS_FVM.py") 文件以生成径向网格并初始化变量：

```bash
python NS_FVM.py
```

2. 运行 [`solver.py`](command:_github.copilot.openRelativePath?%5B%7B%22scheme%22%3A%22file%22%2C%22authority%22%3A%22%22%2C%22path%22%3A%22%2Fd%3A%2FDesktop%2F%E4%B8%AD%E5%9B%BD%E5%BB%BA%E6%9D%90%2FCFD%E7%83%AD%E8%9C%A1%2FInsulating_layer_PISO%2Fsolver.py%22%2C%22query%22%3A%22%22%2C%22fragment%22%3A%22%22%7D%2C%221e544a3e-16c3-4e92-bcb9-b8326c54e279%22%5D "d:\Desktop\中国建材\CFD热蜡\Insulating_layer_PISO\solver.py") 文件以求解 Navier-Stokes 方程并生成结果图像：

```bash
python solver.py
```

3. 生成的图像将保存在项目目录中，包括压力场分布图、速度场分布图和温度分布图等。

## 示例输出

以下是一些示例输出图像：

- [`Pressure.png`](command:_github.copilot.openRelativePath?%5B%7B%22scheme%22%3A%22file%22%2C%22authority%22%3A%22%22%2C%22path%22%3A%22%2Fd%3A%2FDesktop%2F%E4%B8%AD%E5%9B%BD%E5%BB%BA%E6%9D%90%2FCFD%E7%83%AD%E8%9C%A1%2FInsulating_layer_PISO%2FPressure.png%22%2C%22query%22%3A%22%22%2C%22fragment%22%3A%22%22%7D%2C%221e544a3e-16c3-4e92-bcb9-b8326c54e279%22%5D "d:\Desktop\中国建材\CFD热蜡\Insulating_layer_PISO\Pressure.png")：压力场分布图。
- [`RadialMesh.png`](command:_github.copilot.openRelativePath?%5B%7B%22scheme%22%3A%22file%22%2C%22authority%22%3A%22%22%2C%22path%22%3A%22%2Fd%3A%2FDesktop%2F%E4%B8%AD%E5%9B%BD%E5%BB%BA%E6%9D%90%2FCFD%E7%83%AD%E8%9C%A1%2FInsulating_layer_PISO%2FRadialMesh.png%22%2C%22query%22%3A%22%22%2C%22fragment%22%3A%22%22%7D%2C%221e544a3e-16c3-4e92-bcb9-b8326c54e279%22%5D "d:\Desktop\中国建材\CFD热蜡\Insulating_layer_PISO\RadialMesh.png")：径向网格的可视化。
- [`Temperature.png`](command:_github.copilot.openRelativePath?%5B%7B%22scheme%22%3A%22file%22%2C%22authority%22%3A%22%22%2C%22path%22%3A%22%2Fd%3A%2FDesktop%2F%E4%B8%AD%E5%9B%BD%E5%BB%BA%E6%9D%90%2FCFD%E7%83%AD%E8%9C%A1%2FInsulating_layer_PISO%2FTemperature.png%22%2C%22query%22%3A%22%22%2C%22fragment%22%3A%22%22%7D%2C%221e544a3e-16c3-4e92-bcb9-b8326c54e279%22%5D "d:\Desktop\中国建材\CFD热蜡\Insulating_layer_PISO\Temperature.png")：温度分布图。
- [`VelocityQuiver.png`](command:_github.copilot.openRelativePath?%5B%7B%22scheme%22%3A%22file%22%2C%22authority%22%3A%22%22%2C%22path%22%3A%22%2Fd%3A%2FDesktop%2F%E4%B8%AD%E5%9B%BD%E5%BB%BA%E6%9D%90%2FCFD%E7%83%AD%E8%9C%A1%2FInsulating_layer_PISO%2FVelocityQuiver.png%22%2C%22query%22%3A%22%22%2C%22fragment%22%3A%22%22%7D%2C%221e544a3e-16c3-4e92-bcb9-b8326c54e279%22%5D "d:\Desktop\中国建材\CFD热蜡\Insulating_layer_PISO\VelocityQuiver.png")：速度场矢量图。
- [`XVelocity.png`](command:_github.copilot.openRelativePath?%5B%7B%22scheme%22%3A%22file%22%2C%22authority%22%3A%22%22%2C%22path%22%3A%22%2Fd%3A%2FDesktop%2F%E4%B8%AD%E5%9B%BD%E5%BB%BA%E6%9D%90%2FCFD%E7%83%AD%E8%9C%A1%2FInsulating_layer_PISO%2FXVelocity.png%22%2C%22query%22%3A%22%22%2C%22fragment%22%3A%22%22%7D%2C%221e544a3e-16c3-4e92-bcb9-b8326c54e279%22%5D "d:\Desktop\中国建材\CFD热蜡\Insulating_layer_PISO\XVelocity.png")：X方向速度场分布图。
- [`YVelocity.png`](command:_github.copilot.openRelativePath?%5B%7B%22scheme%22%3A%22file%22%2C%22authority%22%3A%22%22%2C%22path%22%3A%22%2Fd%3A%2FDesktop%2F%E4%B8%AD%E5%9B%BD%E5%BB%BA%E6%9D%90%2FCFD%E7%83%AD%E8%9C%A1%2FInsulating_layer_PISO%2FYVelocity.png%22%2C%22query%22%3A%22%22%2C%22fragment%22%3A%22%22%7D%2C%221e544a3e-16c3-4e92-bcb9-b8326c54e279%22%5D "d:\Desktop\中国建材\CFD热蜡\Insulating_layer_PISO\YVelocity.png")：Y方向速度场分布图。
- [`MassResidual.png`](command:_github.copilot.openRelativePath?%5B%7B%22scheme%22%3A%22file%22%2C%22authority%22%3A%22%22%2C%22path%22%3A%22%2Fd%3A%2FDesktop%2F%E4%B8%AD%E5%9B%BD%E5%BB%BA%E6%9D%90%2FCFD%E7%83%AD%E8%9C%A1%2FInsulating_layer_PISO%2FMassResidual.png%22%2C%22query%22%3A%22%22%2C%22fragment%22%3A%22%22%7D%2C%221e544a3e-16c3-4e92-bcb9-b8326c54e279%22%5D "d:\Desktop\中国建材\CFD热蜡\Insulating_layer_PISO\MassResidual.png"): 质量残差图。
- [`XvelocityResidual.png`](command:_github.copilot.openRelativePath?%5B%7B%22scheme%22%3A%22file%22%2C%22authority%22%3A%22%22%2C%22path%22%3A%22%2Fd%3A%2FDesktop%2F%E4%B8%AD%E5%9B%BD%E5%BB%BA%E6%9D%90%2FCFD%E7%83%AD%E8%9C%A1%2FInsulating_layer_PISO%2FXvelocityResidual.png%22%2C%22query%22%3A%22%22%2C%22fragment%22%3A%22%22%7D%2C%221e544a3e-16c3-4e92-bcb9-b8326c54e279%22%5D "d:\Desktop\中国建材\CFD热蜡\Insulating_layer_PISO\XvelocityResidual.png"): X方向速度残差图。
- [`YvelocityResidual.png`](command:_github.copilot.openRelativePath?%5B%7B%22scheme%22%3A%22file%22%2C%22authority%22%3A%22%22%2C%22path%22%3A%22%2Fd%3A%2FDesktop%2F%E4%B8%AD%E5%9B%BD%E5%BB%BA%E6%9D%90%2FCFD%E7%83%AD%E8%9C%A1%2FInsulating_layer_PISO%2FYvelocityResidual.png%22%2C%22query%22%3A%22%22%2C%22fragment%22%3A%22%22%7D%2C%221e544a3e-16c3-4e92-bcb9-b8326c54e279%22%5D "d:\Desktop\中国建材\CFD热蜡\Insulating_layer_PISO\YvelocityResidual.png"): Y方向速度残差图。

## 贡献指南

如果你想为该项目做出贡献，请遵循以下步骤：

1. Fork 该项目。
2. 创建一个新的分支：

```bash
git checkout -b feature-branch
```

3. 提交你的更改：

```bash
git commit -am 'Add some feature'
```

4. 推送到分支：

```bash
git push origin feature-branch
```

5. 创建一个新的 Pull Request。

## 许可证

该项目使用 MIT 许可证。详情请参阅 LICENSE 文件。