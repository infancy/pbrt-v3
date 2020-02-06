#include "..."

// Chapter1
class Shape;
class Material;
class Primitive;
class Accelerator;
class Spectrum;
class Camera;
class Sampler;
class Filter;
class Film;
class BSDF;
class Texture;
class Meida;
class Light;
class Integraotr; // Renderer

// 还有主渲染的代码

// Chapter2
class VectorN{ ... };
class Normal3{ ... };
class Ray{ ... };
class RayDifference{ ... };

// Chapter3
class Shape
{
    bool Interaction();
};

// Chapter4
class Primitive
{
    bool Interaction();
};

// ...

// 最后再回过头来看看主渲染流程