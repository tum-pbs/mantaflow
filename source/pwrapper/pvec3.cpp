/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Vec3 class extension for python
 *
 ******************************************************************************/

#include <math.h>
#include <string>
#include <sstream>
#include <limits>
#include "pythonInclude.h"
#include "vectorbase.h"
#include "structmember.h"
#include "pclass.h"

using namespace std;

namespace Manta {
  
extern PyTypeObject PbVec3Type;

struct PbVec3 {
    PyObject_HEAD
    float data[3];
};

static void PbVec3Dealloc(PbVec3* self) {
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static PyObject * PbVec3New(PyTypeObject *type, PyObject *args, PyObject *kwds) {    
    return type->tp_alloc(type, 0);
}

static int PbVec3Init(PbVec3 *self, PyObject *args, PyObject *kwds) {
    
    float x1 = numeric_limits<float>::quiet_NaN(), x2=x1, x3=x1;
    if (!PyArg_ParseTuple(args,"|fff",&x1, &x2, &x3))
        return -1;
    
    if (!c_isnan(x1)) {
        self->data[0] = x1;
        if (!c_isnan(x2) && !c_isnan(x3)) {
            self->data[1] = x2;
            self->data[2] = x3;
        } else {
            self->data[1] = x1;
            self->data[2] = x1;
        }
    } else {
        self->data[0] = 0;
        self->data[1] = 0;
        self->data[2] = 0;
    }
    return 0;
}

static PyObject* PbVec3Repr(PbVec3* self) {
    Manta::Vec3 v(self->data[0], self->data[1], self->data[2]);
    return PyUnicode_FromFormat(v.toString().c_str());
}

static PyMemberDef PbVec3Members[] = {
    {(char*)"x", T_FLOAT, offsetof(PbVec3, data), 0, (char*)"X"},
    {(char*)"y", T_FLOAT, offsetof(PbVec3, data)+sizeof(float), 0, (char*)"Y"},
    {(char*)"z", T_FLOAT, offsetof(PbVec3, data)+sizeof(float)*2, 0, (char*)"Z"},
    {NULL}  // Sentinel
};

static PyMethodDef PbVec3Methods[] = {
    //{"name", (PyCFunction)Noddy_name, METH_NOARGS, "Return the name, combining the first and last name" },
    {NULL}  // Sentinel
};

// operator overloads
    
inline PyObject* PbNew(const Vec3& a) { 
    PbVec3* obj = (PbVec3*)PbVec3New(&PbVec3Type,0,0); 
    obj->data[0] = a.x;
    obj->data[1] = a.y;
    obj->data[2] = a.z;
    return (PyObject*)obj;
}

#define CONVERTVEC(obj) \
    Vec3 v##obj; \
    if (PyObject_TypeCheck(obj, &PbVec3Type)) \
        v##obj = Vec3(&( ((PbVec3*)obj)->data[0])); \
    else if (PyFloat_Check(obj)) \
        v##obj = Vec3(PyFloat_AsDouble(obj)); \
    else if (PyLong_Check(obj)) \
        v##obj = Vec3(PyLong_AsDouble(obj)); \
    else { \
        Py_INCREF(Py_NotImplemented); \
        return Py_NotImplemented; \
    } \

#define OPHEADER \
    if (!PyObject_TypeCheck(a, &PbVec3Type) && !PyObject_TypeCheck(b, &PbVec3Type)) { \
        Py_INCREF(Py_NotImplemented); \
        return Py_NotImplemented; \
    } \
    CONVERTVEC(a) \
    CONVERTVEC(b)

#define OPHEADER1 \
    if (!PyObject_TypeCheck(a, &PbVec3Type)) { \
        Py_INCREF(Py_NotImplemented); \
        return Py_NotImplemented; \
    } \
    CONVERTVEC(a)
    
PyObject* PbVec3Add(PyObject* a, PyObject* b) {
    OPHEADER
    return PbNew(va+vb);
}

PyObject* PbVec3Sub(PyObject* a, PyObject* b) {
    OPHEADER
    return PbNew(va-vb);
}

PyObject* PbVec3Mult(PyObject* a, PyObject* b) {
    OPHEADER
    return PbNew(va*vb);
}

PyObject* PbVec3Div(PyObject* a, PyObject* b) {
    OPHEADER
    return PbNew(va/vb);
}

PyObject* PbVec3Negative(PyObject* a) {
    OPHEADER1
    return PbNew(-va);
}


static PyNumberMethods PbVec3NumberMethods = {
     (binaryfunc)PbVec3Add,          // binaryfunc nb_add;
     (binaryfunc)PbVec3Sub,          // binaryfunc nb_sub;
     (binaryfunc)PbVec3Mult,         // binaryfunc nb_mult;
     0,          // binaryfunc nb_remainder;
     0,          // binaryfunc nb_divmod;
     0,          // ternaryfunc nb_power;
     (unaryfunc)PbVec3Negative,          // unaryfunc nb_negative;
     0,          // unaryfunc nb_positive;
     0,          // unaryfunc nb_absolute;
     0,          // inquiry nb_bool;
     0,          // unaryfunc nb_invert;
     0,          // binaryfunc nb_lshift;
     0,          // binaryfunc nb_rshift;
     0,          // binaryfunc nb_and;
     0,          // binaryfunc nb_xor;
     0,          // binaryfunc nb_or;
     0,          // unaryfunc nb_int;
     0,          // void *nb_reserved;
     0,          // unaryfunc nb_float;

     0,          // binaryfunc nb_inplace_add;
     0,          // binaryfunc nb_inplace_subtract;
     0,          // binaryfunc nb_inplace_multiply;
     0,          // binaryfunc nb_inplace_remainder;
     0,          // ternaryfunc nb_inplace_power;
     0,          // binaryfunc nb_inplace_lshift;
     0,          // binaryfunc nb_inplace_rshift;
     0,          // binaryfunc nb_inplace_and;
     0,          // binaryfunc nb_inplace_xor;
     0,          // binaryfunc nb_inplace_or;

     0,          // binaryfunc nb_floor_divide;
     (binaryfunc)PbVec3Div,          // binaryfunc nb_true_divide;
     0,          // binaryfunc nb_inplace_floor_divide;
     0,          // binaryfunc nb_inplace_true_divide;

     0           // unaryfunc nb_index;
};

PyTypeObject PbVec3Type = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "manta.vec3",             /* tp_name */
    sizeof(PbVec3),             /* tp_basicsize */
    0,                         /* tp_itemsize */
    (destructor)PbVec3Dealloc, /* tp_dealloc */
    0,                         /* tp_print */
    0,                         /* tp_getattr */
    0,                         /* tp_setattr */
    0,                         /* tp_reserved */
    (reprfunc)PbVec3Repr,      /* tp_repr */
    &PbVec3NumberMethods,      /* tp_as_number */
    0,                         /* tp_as_sequence */
    0,                         /* tp_as_mapping */
    0,                         /* tp_hash  */
    0,                         /* tp_call */
    0,                         /* tp_str */
    0,                         /* tp_getattro */
    0,                         /* tp_setattro */
    0,                         /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,   /* tp_flags */
    "float vector type",           /* tp_doc */
    0,                     /* tp_traverse */
    0,                     /* tp_clear */
    0,                     /* tp_richcompare */
    0,                     /* tp_weaklistoffset */
    0,                     /* tp_iter */
    0,                     /* tp_iternext */
    PbVec3Methods,             /* tp_methods */
    PbVec3Members,             /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)PbVec3Init,      /* tp_init */
    0,                         /* tp_alloc */
    PbVec3New,                 /* tp_new */
};

void PbVecInitialize(PyObject* module) {
    if (PyType_Ready(&PbVec3Type) < 0) throw Manta::Error("can't initialize Vec3 type");
    
    Py_INCREF(&PbVec3Type);
    PyModule_AddObject(module, "vec3", (PyObject *)&PbVec3Type);
}
PB_REGISTER_EXTERNAL(PbVecInitialize);

} // namespace