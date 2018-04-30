//Copyright (c) 2013-2018, The MercuryDPM Developers Team. All rights reserved.
//For the list of developers, see <http://www.MercuryDPM.org/Team>.
//
//Redistribution and use in source and binary forms, with or without
//modification, are permitted provided that the following conditions are met:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name MercuryDPM nor the
//    names of its contributors may be used to endorse or promote products
//    derived from this software without specific prior written permission.
//
//THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
//ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
//WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//DISCLAIMED. IN NO EVENT SHALL THE MERCURYDPM DEVELOPERS TEAM BE LIABLE FOR ANY
//DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
//(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
//LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
//ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
//(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
//SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

/** Implementation of a 3D vector (by Vitaliy).
 */
///Modifications
/// 21:9:2009 - Added the inclusion guards and some include objects
/// \todo Need to generalize this to n-dimensional vectors of any type
#ifndef VECTOR_H
#define VECTOR_H

#include <cmath>
#include <sstream>
#include <iostream>
#include <stdlib.h> 

#include "GeneralDefine.h"
#include "Logger.h"

template <unsigned int N>
class SmallVector;

/*!
 * \class Vec3D
 * \brief
 */
class Vec3D
{
public:
    
    /*!
     * \brief the vector components
     */
    /*
     * \todo: Make these private.
     * These should be private so we can implement things like a cvec etc.
     * Use getters / setters.
     */
    // private:
    Mdouble X, Y, Z;

    /*!
     * \brief constructor
     */
    Vec3D() {setZero();}
    
    Vec3D(const SmallVector<3>& vector);
    
    /*!
     * \brief Alternative constructor, taking the three elements as arguments
     */
    Vec3D(const Mdouble x, const Mdouble y, const Mdouble z) { X = x; Y = y; Z = z; }
    
    /*!
     * \brief Sets all elements to zero
     */
    void setZero();

    /*!
     * \brief Sets all elements to NaN
     */
    void setNaN();

    /*!
     * \brief Checks if ALL elements are zero
     */
    bool isZero() const;

    /*!
     * \brief Checks if ALL elements are zero
     */
    bool isNaN() const;

    /*!
     * \brief Adds another vector
     */
    Vec3D operator+(const Vec3D& a) const;

    /*!
     * \brief Subtracts another vector
     * \details Subtracts vector from itself
     * \param[in] a     vector to be subtracted
     * \return          resulting vector
     */
    Vec3D operator-(const Vec3D& a) const {
        return Vec3D(X - a.X, Y - a.Y, Z - a.Z);
    };

    /*!
     * \brief Multiplies by a scalar
     * \details Multiplies each element with a scalar
     * \param[in] a     the scalar to be multiplied with
     * \return          the resulting vector
     */
    Vec3D operator*(const Mdouble a) const {
        return Vec3D(X * a, Y * a, Z * a);
    }

    /*!
     * \brief Divides by a scalar
     */
    Vec3D operator/(const Mdouble a) const;
    
    /*!
     * \brief Adds another vector
     * \details Adds a vector to itself
     * \param[in] a     vector to be added
     * \return          (reference to) itself, i.e. resulting vector
     */
    Vec3D& operator+=(const Vec3D& a) {
        X += a.X;
        Y += a.Y;
        Z += a.Z;
        return *this;
    }
    
    /*!
     * \brief Subtracts another vector
     */
    Vec3D& operator-=(const Vec3D& a);
    
    /*!
     * \brief Multiplies by a scalar
     */
    Vec3D& operator*=(const Mdouble a);
    
    /*!
     * \brief Divides by a scalar
     */
    Vec3D& operator/=(const Mdouble a);
    
    /*!
     * \brief Calculates the dot product of two Vec3D: \f$ a \cdot b\f$
     */
    static Mdouble dot(const Vec3D& a, const Vec3D& b);
    
    /*!
     * \brief Calculates the pointwise maximum of two Vec3D
     */
    static Vec3D max(const Vec3D& a, const Vec3D& b);
    
    /*!
     * \brief Calculates the pointwise minimum of two Vec3D
     */
    static Vec3D min(const Vec3D& a, const Vec3D& b);
    
    /*!
     * \brief Calculates the pointwise square of a Vec3D
     */
    static Vec3D square(const Vec3D& a);
    
    /*!
     * \brief Makes this Vec3D unit length
     */
    void normalize();
    
    /*!
     * \brief Make this Vec3D a certain length 
     */
    void setLength(Mdouble length);
    
    /*!
     * \brief Calculates the pointwise square root of a Vec3D
     */
    static Vec3D sqrt(const Vec3D& a);
    
    /*!
     * \brief Calculates the cross product of two Vec3D:  \f$ a \times b\f$
     */
    static Vec3D cross(const Vec3D& a, const Vec3D& b);
    
    /*!
     * \brief Calculates the distance between two Vec3D: \f$ \sqrt{\left(a-b\right) \cdot \left(a-b\right)} \f$
     */
    static Mdouble getDistance(const Vec3D& a, const Vec3D& b);
    
    /*!
     * \brief Calculates the squared distance between two Vec3D: \f$ \left(a-b\right) \cdot \left(a-b\right) \f$
     */
    static Mdouble getDistanceSquared(const Vec3D& a, const Vec3D& b);
    
    /*!
     * \brief Calculates the length of a Vec3D: \f$ \sqrt{a\cdot a} \f$
     */
    static Mdouble getLength(const Vec3D& a);

    /*!
     * \brief Calculates the squared length of a Vec3D: \f$ a\cdot a \f$
     * \details Calculates the square of the length of a given vector.
     * NB: this is a STATIC function!
     * \param[in] a     the vector.
     * \return          the square of the length of the argument.
     */
    static Mdouble getLengthSquared(const Vec3D& a) {
        return (a.X * a.X + a.Y * a.Y + a.Z * a.Z);
    }
    
    /*!
     * \brief Calculates the length of this Vec3D: \f$ \sqrt{a\cdot a} \f$
     */
    Mdouble getLength() const;
    
    /*!
     * \brief Calculates the squared length of this Vec3D: \f$ a\cdot a \f$
     */
    Mdouble getLengthSquared() const;
    
    /*!
     * \brief Returns the requested component of this Vec3D
     */
    Mdouble getComponent(const int index) const;
    
    /*!
     * \brief Sets the requested component of this Vec3D to the requested value
     */
    void setComponent(const int index, const double val);
    
    /*!
     * \brief RW reference to X
     */
    inline Mdouble& x() { return X; }
    
    /*!
     * \brief RO reference to X
     */
    inline Mdouble x() const { return X; }
    
    /*!
     * \brief RW reference to Y
     */
    inline Mdouble& y() { return Y; }
    
    /*!
     * \brief RO reference to Y
     */
    inline Mdouble y() const { return Y; }
    
    /*!
     * \brief RW reference to Z
     */
    inline Mdouble& z() { return Z; }
    
    /*!
     * \brief RO reference to Z
     */
    inline Mdouble z() const { return Z; }
    
    inline void setX(Mdouble x) { X = x; }
    inline void setY(Mdouble y) { Y = y; }
    inline void setZ(Mdouble z) { Z = z; }
    inline Mdouble getX() { return X; }
    inline Mdouble getY() { return Y; }
    inline Mdouble getZ() { return Z; }
    inline void set(Mdouble x, Mdouble y, Mdouble z) { X = x; Y = y; Z = z; }

    /*!
     * \brief Returns the square of the radial cylindrical coordinate, r^2=x^2+y^2.
     */
    Mdouble getRadialCoordinateSquared() const;

    /*!
     * \brief Returns the square of the radial cylindrical coordinate, r=sqrt(x^2+y^2).
     */
    Mdouble getRadialCoordinate() const;

    /*!
     * \brief Returns the representation of this Vec3D in cylindrical coordinates
     */
    Vec3D getCylindricalCoordinates() const;

    /*!
     * \brief Returns the representation of this Vec3D in cylindrical coordinates
     */
    Vec3D getFromCylindricalCoordinates() const;

    /*!
     * \brief Returns this vector field at point p to cylindrical coordinates
     */
    Vec3D getCylindricalTensorField(const Vec3D& position) const;

    /*!
     * \brief Checks if the length this Vec3D is equal the length of other with a certain tolerance
     */
    bool isEqualTo(const Vec3D& other, const double tol) const;

    /*!
     * \brief Returns a unit Vec3D based on a.
     */
    static Vec3D getUnitVector(const Vec3D& a);
    
    /*!
     * \brief Adds elements to an output stream
     */
    friend std::ostream& operator<<(std::ostream& os, const Vec3D& a);
    
    /*!
     * \brief Adds elements to an input stream
     */
    friend std::istream& operator>>(std::istream& is, Vec3D& a);

    /*!
     * \brief Subtracts a vector 
     */
    friend Vec3D operator-(const Vec3D& a);
    
    /*!
     * \brief Multiplies all elements by a scalar
     */
    friend Vec3D operator*(const Mdouble a, const Vec3D& b);    
    
    
};

#endif
