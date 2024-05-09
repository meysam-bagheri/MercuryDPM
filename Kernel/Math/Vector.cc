//Copyright (c) 2013-2023, The MercuryDPM Developers Team. All rights reserved.
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

#include "Vector.h"
#include "SmallVector.h"
#include "ExtendedMath.h"

/*!
 * \details Alternative constructor, that constructs a Vec3D from a SmallVector size 3
 * \param[in] vector Small vector that should be copied
 */
Vec3D::Vec3D(const SmallVector<3>& vector)
{
    X = vector[0];
    Y = vector[1];
    Z = vector[2];
}

/*!
 * \details Sets each element to zero.
 */
void Vec3D::setZero()
{
    X = 0.0;
    Y = 0.0;
    Z = 0.0;
}

/*!
 * \details Sets each element to zero.
 */
void Vec3D::setNaN()
{
    X = constants::NaN;
    Y = constants::NaN;
    Z = constants::NaN;
}

/*!
 * \details Checks if ALL elements are zero
 * \return          TRUE if ALL elements are zero
 */
bool Vec3D::isNaN() const
{
    return std::isnan(X) || std::isnan(Y) || std::isnan(Z);
}

/*!
 * \details Calculates the dot product of two vectors.
 * NB: this is a STATIC function!
 * \param[in] a     the first vector
 * \param[in] b     the second vector 
 * \return          the resulting scalar
 */
Mdouble Vec3D::dot(const Vec3D& a, const Vec3D& b)
{
    return a.X * b.X + a.Y * b.Y + a.Z * b.Z;
}

/*!
 * \details Calculates the pointwise maximum of two vectors.
 * NB: this is a STATIC function!
 * \param[in] a     the first vector
 * \param[in] b     the second vector 
 * \return          The resulting vector, in which each element is the maximum 
 *                  of the equivalent elements of the arguments
 */
Vec3D Vec3D::max(const Vec3D& a, const Vec3D& b)
{
    return Vec3D(std::max(a.X, b.X), std::max(a.Y, b.Y), std::max(a.Z, b.Z));
}

/*!
 * \details Calculates the pointwise minimum of two vectors.
 * NB: this is a STATIC function!
 * \param[in] a     the first vector
 * \param[in] b     the second vector 
 * \return          The resulting vector, in which each element is the minimum 
 *                  of the equivalent elements of the arguments
 */
Vec3D Vec3D::min(const Vec3D& a, const Vec3D& b)
{
    return Vec3D(std::min(a.X, b.X), std::min(a.Y, b.Y), std::min(a.Z, b.Z));
}

/*!
 * \details Calculates the pointwise square of the vector.
 * NB: this is a STATIC function!
 * \param[in] a     the vector to be squared.
 * \return          the resulting vector, of which each element is the square of
 *                  the equivalent element of the argument.
 */
Vec3D Vec3D::square(const Vec3D& a)
{
    return Vec3D(a.X * a.X, a.Y * a.Y, a.Z * a.Z);
}

/*!
 * \details Normalises the vector, i.e. divides all elements by the vectors length
 * (resulting in a vector in the same direction, but with unit length).
 */
void Vec3D::normalise()
{
    Mdouble length2 = this->getLengthSquared();
    if (length2 == 0)
    {
        logger(ERROR, "Normalizing a vector of length 0");
    }
    *this /= std::sqrt(length2);
}

/*!
 * \details Sets the length of the vector to a given scalar (while maintaining the
 * direction).
 * \param[in] length    the length to be set
 */
void Vec3D::setLength(Mdouble length)
{
    this->normalise();
    *this *= length;
}

/*!
 * \details Calculates the pointwise square root of a given vector.
 * NB: this is a STATIC function!
 * \param[in] a     the vector to be pointwise square rooted
 * \return          the resulting vector, of which each element is the square root
 *                  of the equivalent element of the argument.
 */
Vec3D Vec3D::sqrt(const Vec3D& a)
{
    return Vec3D(std::sqrt(a.X), std::sqrt(a.Y), std::sqrt(a.Z));
}

/*!
 * \details Calculates the cross product of two vectors
 * NB: this is a STATIC function!
 * \param[in] a     the first vector
 * \param[in] b     the second vector 
 * \return          the cross product of the arguments
 */
Vec3D Vec3D::cross(const Vec3D& a, const Vec3D& b)
{
    return Vec3D(a.Y * b.Z - a.Z * b.Y, a.Z * b.X - a.X * b.Z, a.X * b.Y - a.Y * b.X);
}

/*!
 * \details Calculates the distance (i.e. the length of the difference) between two vectors
 * NB: this is a STATIC function!
 * \param[in] a     the first vector
 * \param[in] b     the second vector 
 * \return          the distance between the two arguments.
 */
Mdouble Vec3D::getDistance(const Vec3D& a, const Vec3D& b)
{
    return std::sqrt(getDistanceSquared(a, b));
}

/*!
 * \details Calculates the square of the length of itself
 * \return              the square of the length of this vector
 */
Mdouble Vec3D::getLengthSquared() const
{
    return (X * X + Y * Y + Z * Z);
}

/*!
 * \details returns the vector element belonging to the given index.
 * \param[in] index     the index of interest (should be 0, 1 or 2)
 * \return              the value of the vector element belonging to the given index
 */
Mdouble Vec3D::getComponent(const int index) const
{
    switch (index)
    {
        case 0:
            return X;
        case 1:
            return Y;
        case 2:
            return Z;
        default:
            logger(ERROR, "[Vector::getComponent] Index = %, which is too high for a 3D vector (should be 0-2).",
                   index);
            return 0;
    }
}

/*!
 * \details Sets the element of the vector belonging to the first argument 
 * (index) to the value given in the second argument (val).
 * \param[in] index     index of element of interest, 
 * \param[in] val       value to be set
 */
void Vec3D::setComponent(const int index, const double val)
{
    switch (index)
    {
        case 0:
            X = val;
            break;
        case 1:
            Y = val;
            break;
        case 2:
            Z = val;
            break;
        default:
            logger(ERROR, "[Vector::setComponent] Index = %, which is too high for a 3D vector (should be 0-2).",
                   index);
    }
}


Mdouble Vec3D::getRadialCoordinateSquared() const
{
    return X * X + Y * Y;
}

Mdouble Vec3D::getRadialCoordinate() const
{
    return std::sqrt(X * X + Y * Y);
}

/*!
 * \details Transforms the (Cartesian) vector to spherical coordinates
 * \return              Transformed vector
 */
Vec3D Vec3D::getSphericalCoordinates() const
{
    Vec3D sphericalCoord;

    // Calculate radial distance
    sphericalCoord.X = std::sqrt(X * X + Y * Y + Z * Z);

    // Calculate azimuthal angle (theta)
    sphericalCoord.Y = std::acos(Z / sphericalCoord.X);

    // Calculate polar angle (phi)
    if (X == 0.0)
    {
        if (Y >= 0.0)
        {
            sphericalCoord.Z = constants::pi / 2.0;
        }
        else
        {
            sphericalCoord.Z = -constants::pi / 2.0;
        }
    }
    else
    {
        sphericalCoord.Z = std::atan2(Y, X);
        // Adjust the range to [0, 2*pi)
        if (sphericalCoord.Z < 0.0)
        {
            sphericalCoord.Z += 2.0 * constants::pi;
        }
    }

    return sphericalCoord;
}

/*!
 * \details Transforms the (Cartesian) vector to cylindrical coordinates
 * \return              Transformed vector
 */
Vec3D Vec3D::getCylindricalCoordinates() const
{
    return Vec3D(std::sqrt(X * X + Y * Y), std::atan2(Y, X), Z);
}

/*!
 * \details Transforms the (Cartesian) vector to cylindrical coordinates.
 * See https://en.wikipedia.org/wiki/Vector_fields_in_cylindrical_and_spherical_coordinates
 * \return              Transformed vector
 */
Vec3D Vec3D::getCylindricalTensorField(const Vec3D& p) const
{
    //define sin(A)=y/r, cos(A)=x/r
    Mdouble r = std::sqrt(p.X * p.X + p.Y * p.Y);
    Mdouble s = p.Y / r;
    Mdouble c = p.X / r;
    if (r == 0)
    {
        s = 0;
        c = 1;
    }
    return Vec3D(X * c + Y * s, -X * s + Y * c, Z);
}

/*!
 * \details Transforms the (cylindrical) vector to cartesian coordinates
 * \return              Transformed vector
 */
Vec3D Vec3D::getFromCylindricalCoordinates() const
{
    ///\todo
    return Vec3D(X * std::cos(Y), X * std::sin(Y), Z);
}

/*!
 * \details Checks if the length of the vector is equal to the one given in the 
 * first argument (other), with a tolerance given in the second argument (tol).
 * \param[in] other     the 3D vector to check against
 * \param[in] tol       the tolerance
 * \return              returns TRUE if the difference between the lengths of this 
 *                      vector and that given in the first argument (other) is smaller than the 
 *                      given tolerance.
 */
bool Vec3D::isEqualTo(const Vec3D& other, const double tol) const
{
    if ((Vec3D::getLengthSquared(*this - other)) <= tol * tol)
    {
        return true;
    }
    else
    {
        return false;
    }
}

// 	void ConvertToCylindricalCoordinates()
// 	{
// 		double R = sqrt(X*X+Y*Y); Y = atan2(Y,X); X = R; return;
// 	}
// 
// 	void ConvertFromCylindricalCoordinates()
// 	{
// 		double Xnew = X * cos(Y); Y = X * sin(Y); X = Xnew; return;
// 	}

/*!
 * \details Calculates the length of this vector
 * \return          the (scalar) length of this vector
 */
Mdouble Vec3D::getLength() const
{
    return std::sqrt(getLengthSquared());
}

/*!
 * \details Calculates the length of a given vector
 * NB: this is a STATIC function!
 * \param[in] a     vector to be measured.
 * \return          length of the argument.
 */
Mdouble Vec3D::getLength(const Vec3D& a)
{
    return a.getLength();
}

/*!
 * \details Calculates the unit vector of a given vector (unless it is a vector 
 * with zero length; in that case it returns a 3D vector with each element equal
 * to zero).
 * NB: this is a STATIC function!
 * \param[in] a     the vector of interest
 * \return          unit vector in the direction of the argument (unless the 
 *                  argument has length zero; in that case a zero-vector).
 */
Vec3D Vec3D::getUnitVector(const Vec3D& a)
{
    Mdouble Length2 = a.getLengthSquared();
    if (Length2 != 0.0)
        return a / std::sqrt(Length2);
    else
        return Vec3D(0, 0, 0);
}

/*!
 * \details Adds all elements of the vector to an output stream.
 * NB: this is a global function and a friend of the Vec3D class!
 * \param[in] os    the output stream, 
 * \param[in] a     The vector of interest
 * \return          the output stream with vector elements added
 */
std::ostream& operator<<(std::ostream& os, const Vec3D& a)
{
    os << a.X << ' ' << a.Y << ' ' << a.Z;
    return os;
}

/*!
 * \details Reads all elements of a given vector from an input stream.
 * NB: this is a global function and a friend of the Vec3D class!
 * \param[in,out] is    the input stream
 * \param[in,out] a     the vector to be read in
 * \return              the input stream from which the vector elements were read
 */
std::istream& operator>>(std::istream& is, Vec3D& a)
{
    //TW: clearing the stream avoids the nasty problem that the failbit is set to true if numbers below DBL_MIN=1e-308 are read.
    is >> a.X;
    if (is.fail()) {
        logger(VERBOSE,"Failed to read x-value %", a.X);
        //a.X = 0;
        is.clear();
    }
    is >> a.Y;
    if (is.fail()) {
        logger(VERBOSE,"Failed to read y-value %", a.Y);
        //a.Y = 0;
        is.clear();
    }
    is >> a.Z;
    if (is.fail()) {
        logger(VERBOSE,"Failed to read z-value %", a.Z);
        //a.Z = 0;
        is.clear();
    }
    return is;
}