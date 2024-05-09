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
#include "Math/Quaternion.h"
#include "ExtendedMath.h"

/*!
 * \details Default constructor
 */
Quaternion::Quaternion()
{
    setUnity();
}

/*!
 * \details Alternative constructor, lets you define all four elements.
 * \param[in] q0     the q0-component
 * \param[in] q1     the q1-component
 * \param[in] q2     the q2-component
 * \param[in] q3     the q3-component
 */
Quaternion::Quaternion(const Mdouble q0, const Mdouble q1, const Mdouble q2, const Mdouble q3)
{
    this->q0 = q0;
    this->q1 = q1;
    this->q2 = q2;
    this->q3 = q3;
}

/*!
 * \details Sets q0 to 1, and all other elements to zero.
 */
void Quaternion::setUnity()
{
    q0 = 1.0;
    q1 = 0.0;
    q2 = 0.0;
    q3 = 0.0;
}

/*!
 * \details Adds quaternion to itself
 * \param[in] a     quaternion to be added
 * \return          resulting 3D quaternion
 */
Quaternion Quaternion::operator+(const Quaternion& a) const
{
    return Quaternion(q0 + a.q0, q1 + a.q1, q2 + a.q2, q3 + a.q3);
}

/*!
 * \details Subtracts quaternion from itself
 * \param[in] a     quaternion to be subtracted
 * \return          resulting quaternion
 */
Quaternion Quaternion::operator-(const Quaternion& a) const
{
    return Quaternion(q0 - a.q0, q1 - a.q1, q2 - a.q2, q3 - a.q3);
}

/*!
 * \details Multiplies each element with a scalar
 * \param[in] a     the scalar to be multiplied with
 * \return          the resulting quaternion
 */
Quaternion Quaternion::operator*(const Mdouble a) const
{
    return Quaternion(q0 * a, q1 * a, q2 * a, q3 * a);
}

/*!
 * \details Divides each element by a scalar
 * \param[in] a     the scalar to be divided by
 * \return          resulting quaternion
 */
Quaternion Quaternion::operator/(const Mdouble a) const
{
    return Quaternion(q0 / a, q1 / a, q2 / a, q3 / a);
}

/*!
 * \details Adds a quaternion to itself
 * \param[in] a     quaternion to be added
 * \return          (reference to) itself, i.e. resulting quaternion
 */
Quaternion& Quaternion::operator+=(const Quaternion& a)
{
    q0 += a.q0;
    q1 += a.q1;
    q2 += a.q2;
    q3 += a.q3;
    return *this;
}

/*!
 * \details Subtracts a quaternion from itself
 * \param[in] a     quaternion to be subtracted
 * \return          (reference to) itself, i.e. resulting quaternion
 */
Quaternion& Quaternion::operator-=(const Quaternion& a)
{
    q0 -= a.q0;
    q1 -= a.q1;
    q2 -= a.q2;
    q3 -= a.q3;
    return *this;
}

/*!
 * \details Multiplies each element by a scalar
 * \param[in] a     scalar to be multiplied by
 * \return          (reference to) itself, i.e. resulting quaternion
 */
Quaternion& Quaternion::operator*=(const Mdouble a)
{
    q0 *= a;
    q1 *= a;
    q2 *= a;
    q3 *= a;
    return *this;
}

/*!
 * \details Divides each element by a scalar
 * \param[in] a     scalar to be divided by
 * \return          (reference to) itself, i.e. resulting quaternion
 */
Quaternion& Quaternion::operator/=(const Mdouble a)
{
    q0 /= a;
    q1 /= a;
    q2 /= a;
    q3 /= a;
    return *this;
}

/*!
 * \details Normalises the quaternion, i.e. divides all elements by the quaternions length
 * (resulting in a quaternion in the same direction, but with unit length).
 */
void Quaternion::normalise()
{
    const Mdouble length2 = getLengthSquared();
    if (length2 == 0)
    {
        logger(ERROR, "Normalizing a quaternion of length 0");
    }
    *this /= sqrt(length2);
}

/*!
 * \details Sets the length of the quaternion to a given scalar (while maintaining the
 * direction).
 * \param[in] length    the length to be set
 */
void Quaternion::setLength(Mdouble length)
{
    *this /= this->getLength() * length;
}

/*!
 * \details Calculates the distance (i.e. the length of the difference) between two quaternions
 * NB: this is a STATIC function!
 * \param[in] a     the first quaternion
 * \param[in] b     the second quaternion 
 * \return          the distance between the two arguments.
 */
Mdouble Quaternion::getDistance(const Quaternion& a, const Quaternion& b)
{
    return std::sqrt(getDistanceSquared(a, b));
}

/*!
 * \details Calculates the square of the distance (i.e. the length of the difference)
 * between two quaternions.
 * NB: this is a STATIC function!
 * \param[in] a     the first quaternion
 * \param[in] b     the second quaternion 
 * \return          the square of the distance between the two arguments.
 */
Mdouble Quaternion::getDistanceSquared(const Quaternion& a, const Quaternion& b)
{
    return ((a.q0 - b.q0) * (a.q0 - b.q0) + (a.q1 - b.q1) * (a.q1 - b.q1) + (a.q2 - b.q2) * (a.q2 - b.q2) +
            (a.q3 - b.q3) * (a.q3 - b.q3));
}

/*!
 * \details Calculates the square of the length of a given quaternion.
 * NB: this is a STATIC function!
 * \param[in] a     the quaternion.
 * \return          the square of the length of the argument.
 */
Mdouble Quaternion::getLengthSquared(const Quaternion& a)
{
    return (a.q0 * a.q0 + a.q1 * a.q1 + a.q2 * a.q2 + a.q3 * a.q3);
}

/*!
 * \details Calculates the square of the length of itself
 * \return              the square of the length of this quaternion
 */
Mdouble Quaternion::getLengthSquared() const
{
    return (q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3);
}

/*!
 * \details returns the quaternion element belonging to the given index.
 * \param[in] index     the index of interest (should be 0, 1 or 2)
 * \return              the value of the quaternion element belonging to the given index
 */
Mdouble Quaternion::getComponent(const int index) const
{
    switch (index)
    {
        case 0:
            return q0;
        case 1:
            return q1;
        case 2:
            return q2;
        case 3:
            return q3;
        default:
            logger(ERROR,
                   "[Quaternion::getComponent] Index = %, which is too high for a 3D quaternion (should be 0-2).",
                   index);
            return 0;
    }
}

/*!
 * \details Sets the element of the quaternion belonging to the first argument 
 * (index) to the value given in the second argument (val).
 * \param[in] index     index of element of interest, 
 * \param[in] val       value to be set
 */
void Quaternion::setComponent(const int index, const double val)
{
    switch (index)
    {
        case 0:
            q0 = val;
            break;
        case 1:
            q1 = val;
            break;
        case 2:
            q2 = val;
            break;
        case 3:
            q3 = val;
            break;
        default:
            logger(ERROR,
                   "[Quaternion::setComponent] Index = %, which is too high for a 3D quaternion (should be 0-2).",
                   index);
    }
}

/*!
 * \details Checks if the length of the quaternion is equal to the one given in the 
 * first argument (other), with a tolerance given in the second argument (tol).
 * \param[in] other     the 3D quaternion to check against
 * \param[in] tol       the tolerance
 * \return              returns TRUE if the difference between the lengths of this 
 *                      quaternion and that given in the first argument (other) is smaller than the 
 *                      given tolerance.
 */
bool Quaternion::isEqualTo(const Quaternion& other, const double tol) const
{
    return (getDistance(*this, other) <= tol * tol);
}

/*!
 * \details Calculates the length of this quaternion
 * \return          the (scalar) length of this quaternion
 */
Mdouble Quaternion::getLength() const
{
    return std::sqrt(getLengthSquared());
}

/*!
 * \details Calculates the length of a given quaternion
 * NB: this is a STATIC function!
 * \param[in] a     quaternion to be measured.
 * \return          length of the argument.
 */
Mdouble Quaternion::getLength(const Quaternion& a)
{
    return a.getLength();
}

/*!
 * \details Calculates the unit quaternion of a given quaternion (unless it is a quaternion 
 * with zero length; in that case it returns a 3D quaternion with each element equal
 * to zero).
 * NB: this is a STATIC function!
 * \param[in] a     the quaternion of interest
 * \return          unit quaternion in the direction of the argument (unless the 
 *                  argument has length zero; in that case a zero-quaternion).
 */
Quaternion Quaternion::getUnitQuaternion(const Quaternion& a)
{
    Mdouble Length2 = a.getLengthSquared();
    if (Length2 != 0.0)
        return a / std::sqrt(Length2);
    else
        return Quaternion(1, 0, 0, 0);
}

/*!
 * \details Adds all elements of the quaternion to an output stream.
 * NB: this is a global function and a friend of the Quaternion class!
 * \param[in] os    the output stream, 
 * \param[in] a     The quaternion of interest
 * \return          the output stream with quaternion elements added
 */
std::ostream& operator<<(std::ostream& os, const Quaternion& a)
{
    os << a.q0 << ' ' << a.q1 << ' ' << a.q2 << ' ' << a.q3;
    return os;
}

/*!
 * \details Reads all elements of a given quaternion from an input stream.
 * NB: this is a global function and a friend of the Quaternion class!
 * \param[in,out] is    the input stream
 * \param[in,out] a     the quaternion to be read in
 * \return              the input stream from which the quaternion elements were read
 */
std::istream& operator>>(std::istream& is, Quaternion& a)
{
    is >> a.q0 >> a.q1 >> a.q2 >> a.q3;
    return is;
}

/*!
 * \details Adds a scalar to the elements of given quaternion
 * NB this is a global function and a friend of the Quaternion class. Gets called when
 * addition operation of the form (Mdouble) + (Quaternion) is performed.
 * \param[in] a     the scalar to be added
 * \param[in] b     the quaternion the scalar gets added to.
 * \return          the resulting quaternion.
 */
Quaternion operator+(const Mdouble a, const Quaternion& b)
{
    return Quaternion(b.q0 + a, b.q1 + a, b.q2 + a, b.q3 + a);
}

/*!
 * \details Subtracts each element of a given quaternion from a scalar 
 * NB this is a global function and a friend of the Quaternion class. Gets called when
 * subtraction operation of the form (Mdouble) - (Quaternion) is performed.
 * \param[in] a     the scalar
 * \param[in] b     the quaternion to be subtracted the scalar gets subtracted from.
 * \return          the resulting quaternion.
 */
Quaternion operator-(const Mdouble a, const Quaternion& b)
{
    return Quaternion(a - b.q0, a - b.q1, a - b.q2, a - b.q3);
}

/*!
 * \details Returns the negative of a given quaternion.
 * NB: this is a global function and a friend of the Quaternion class. Gets called when
 * a negation operation of the form - (Quaternion) is performed. 
 * \param[in] a     the quaternion to be negated
 * \return          the negated quaternion
 */
Quaternion operator-(const Quaternion& a)
{
    return Quaternion(-a.q0, -a.q1, -a.q2, -a.q3);
}

/*!
 * \details Multiplies each element of a given quaternion (b) by a given scalar (a).
 * NB: this is a global function and a friend of the Quaternion class. Gets called when
 * a scalar multiplication of the form (Mdouble) * (Quaternion) is performed.
 * \param[in] a     the scalar 
 * \param[in] b     the quaternion
 * \return          the resulting quaternion
 */
Quaternion operator*(const Mdouble a, const Quaternion& b)
{
    return Quaternion(b.q0 * a, b.q1 * a, b.q2 * a, b.q3 * a);
}

///\todo rename to angularVelocityBodyFixedFrameToAngularDisplacement?
Quaternion Quaternion::angularVelocityBodyFixedFrameToAngularDisplacement(Vec3D v) const
{
    return Quaternion(
            -q1 * v.X - q2 * v.Y - q3 * v.Z,
            q0 * v.X - q3 * v.Y + q2 * v.Z,
            q3 * v.X + q0 * v.Y - q1 * v.Z,
            -q2 * v.X + q1 * v.Y + q0 * v.Z);
}

///Given v = &omega * dt, with omega the angular velocity, this computes the change in angular displacement to be added
///in the time integration. This is equivalent to applying the matrix &tilde{C}
Quaternion Quaternion::angularDisplacementTimeDerivative(Vec3D v) const
{
    return 0.5 * Quaternion(
            -q1 * v.X - q2 * v.Y - q3 * v.Z,
            q0 * v.X + q3 * v.Y - q2 * v.Z,
            -q3 * v.X + q0 * v.Y + q1 * v.Z,
            q2 * v.X - q1 * v.Y + q0 * v.Z);
}

void Quaternion::updateAngularDisplacement(Vec3D angularVelocityDt)
{
    *this += angularDisplacementTimeDerivative(angularVelocityDt);
    //const Quaternion q = *this;
    normalise();
    //logger(INFO,"%|%",q,*this);
}


///\todo rename to angularDisplacementToAngularVelocity?
Vec3D Quaternion::applyCInverse(Quaternion q) const
{
    return 2.0 * Vec3D(
            -q1 * q.q0 + q0 * q.q1 - q3 * q.q2 + q2 * q.q3,
            -q2 * q.q0 + q3 * q.q1 + q0 * q.q2 - q1 * q.q3,
            -q3 * q.q0 - q2 * q.q1 + q1 * q.q2 + q0 * q.q3);
}

///\details Get the Euler angles of this quaternion.
/// Example of visualising Euler angles can be found
/// <a href="http://www.euclideanspace.com/maths/geometry/rotations/euler/examples/"> here </a>
Vec3D Quaternion::getEuler() const
{
    Mdouble sinp = 2 * (q0 * q2 - q3 * q1);
    Mdouble pitch;

    if ((std::abs(sinp) < 1))
    {
        pitch = std::asin(sinp);
    }
    else
    {
        pitch = copysign(constants::pi, sinp);
    }

    return Vec3D(
            std::atan2(2 * (q0 * q1 + q2 * q3), 1 - 2 * (q1 * q1 + q2 * q2)),
            pitch,
            std::atan2(2 * (q0 * q3 + q1 * q2), 1 - 2 * (q2 * q2 + q3 * q3)));
}

void Quaternion::setEuler(const Vec3D& e)
{
    Mdouble cp = std::cos(0.5 * e.Y);
    Mdouble sp = std::sin(0.5 * e.Y);
    Mdouble cr = std::cos(0.5 * e.X);
    Mdouble sr = std::sin(0.5 * e.X);
    Mdouble cy = std::cos(0.5 * e.Z);
    Mdouble sy = std::sin(0.5 * e.Z);
    q0 = cr * cp * cy + sr * sp * sy;
    q1 = sr * cp * cy - cr * sp * sy;
    q2 = cr * sp * cy + sr * cp * sy;
    q3 = cr * cp * sy - sr * sp * cy;
}

Mdouble Quaternion::getAngleZ() const
{
    return -std::atan2(2 * (q0 * q3 + q1 * q2), 1 - 2 * (q2 * q2 + q3 * q3));
}

void Quaternion::setAngleZ(Mdouble psi)
{
    //assuming theta=phi=0
    q0 = std::cos(-0.5 * psi);
    q1 = 0;
    q2 = 0;
    q3 = std::sin(-0.5 * psi);
}

Vec3D Quaternion::getAxis() const
{
    const Mdouble q00 = q0*q0;
    const Mdouble q11 = q1*q1;
    const Mdouble q22 = q2*q2;
    const Mdouble q33 = q3*q3;
    const Mdouble q12 = q1*q2;
    const Mdouble q30 = q3*q0;
    const Mdouble q13 = q1*q3;
    const Mdouble q20 = q2*q0;
    return Vec3D(q00 - q33 + q11 - q22, 2.0 * (q12 + q30), 2.0 * (q13 - q20));
}

//retrieves the rotation matrix, often called A in literature.
void Quaternion::getRotationMatrix(SmallMatrix<3, 3>& A) const
{
    
    Mdouble q00 = q0 * q0;
    Mdouble q01 = 2 * q0 * q1;
    Mdouble q02 = 2 * q0 * q2;
    Mdouble q03 = 2 * q0 * q3;
    Mdouble q11 = q1 * q1;
    Mdouble q12 = 2 * q1 * q2;
    Mdouble q13 = 2 * q1 * q3;
    Mdouble q22 = q2 * q2;
    Mdouble q23 = 2 * q2 * q3;
    Mdouble q33 = q3 * q3;
    A(0, 0) = q00 + q11 - q22 - q33;
    A(1, 0) = q12 + q03;
    A(2, 0) = q13 - q02;
    A(0, 1) = q12 - q03;
    A(1, 1) = q00 - q11 + q22 - q33;
    A(2, 1) = q23 + q01;
    A(0, 2) = q13 + q02;
    A(1, 2) = q23 - q01;
    A(2, 2) = q00 - q11 - q22 + q33;
}

/**
 * Defines *one possible* orientation that rotates the x-axis into the direction given by normal
 * \param[in] normal the vector that the x-axis is rotated into.
 */
void Quaternion::setOrientationViaNormal(Vec3D normal)
{
    //if the normal vector cannot be read properly
//    if (normal.getLengthSquared() < 1e-20)
//    {
//        setUnity();
//        return;
//    }
    
    normal.normalise();
    
    if (normal.X <= -1 + 1e-14)
    {
        *this = Quaternion(0, 0, 1, 0);
        return;
    }
    
    Vec3D half = Vec3D(1, 0, 0) + normal;
    q0 = half.X;
    q1 = 0;
    q2 = -half.Z;
    q3 = half.Y;
    normalise();
    //note, it can technically happen that normalising a normalised vector slightly changes the result.
}

//Euler rodriguez
void Quaternion::rotate(Vec3D& position) const
{
    Mdouble q00 = q0 * q0;
    Mdouble q01 = 2 * q0 * q1;
    Mdouble q02 = 2 * q0 * q2;
    Mdouble q03 = 2 * q0 * q3;
    Mdouble q11 = q1 * q1;
    Mdouble q12 = 2 * q1 * q2;
    Mdouble q13 = 2 * q1 * q3;
    Mdouble q22 = q2 * q2;
    Mdouble q23 = 2 * q2 * q3;
    Mdouble q33 = q3 * q3;
    position = Matrix3D(
            q00 + q11 - q22 - q33, q12 - q03, q13 + q02,
            q12 + q03, q00 - q11 + q22 - q33, q23 - q01,
            q13 - q02, q23 + q01, q00 - q11 - q22 + q33) * position;
}

Vec3D Quaternion::rotate(const Vec3D& position) const
{
    Mdouble q00 = q0 * q0;
    Mdouble q01 = 2 * q0 * q1;
    Mdouble q02 = 2 * q0 * q2;
    Mdouble q03 = 2 * q0 * q3;
    Mdouble q11 = q1 * q1;
    Mdouble q12 = 2 * q1 * q2;
    Mdouble q13 = 2 * q1 * q3;
    Mdouble q22 = q2 * q2;
    Mdouble q23 = 2 * q2 * q3;
    Mdouble q33 = q3 * q3;
    return Matrix3D(
            q00 + q11 - q22 - q33, q12 - q03, q13 + q02,
            q12 + q03, q00 - q11 + q22 - q33, q23 - q01,
            q13 - q02, q23 + q01, q00 - q11 - q22 + q33) * position;
}

///Rotate the given vector from the body-fixed angles to the lab-fixed angles. This is the same as multiplying with the
///rotation matrix, A.
void Quaternion::rotate(SmallVector<3>& position) const
{
    SmallMatrix<3, 3> A;
    getRotationMatrix(A);
    position = A * position;
}

///Rotate the given vector from the lab-fixed angles to the body-fixed angles. This is the same as multiplying with the
///inverse/transpose of the rotation matrix, A^T = A^{-1}.
void Quaternion::rotateBack(Vec3D& position) const
{
    Mdouble q00 = q0 * q0;
    Mdouble q01 = 2 * q0 * q1;
    Mdouble q02 = 2 * q0 * q2;
    Mdouble q03 = 2 * q0 * q3;
    Mdouble q11 = q1 * q1;
    Mdouble q12 = 2 * q1 * q2;
    Mdouble q13 = 2 * q1 * q3;
    Mdouble q22 = q2 * q2;
    Mdouble q23 = 2 * q2 * q3;
    Mdouble q33 = q3 * q3;
    position = Matrix3D(
            q00 + q11 - q22 - q33, q12 + q03, q13 - q02,
            q12 - q03, q00 - q11 + q22 - q33, q23 + q01,
            q13 + q02, q23 - q01, q00 - q11 - q22 + q33) * position;
}

Vec3D Quaternion::rotateBack(const Vec3D& position) const
{
    Mdouble q00 = q0 * q0;
    Mdouble q01 = 2 * q0 * q1;
    Mdouble q02 = 2 * q0 * q2;
    Mdouble q03 = 2 * q0 * q3;
    Mdouble q11 = q1 * q1;
    Mdouble q12 = 2 * q1 * q2;
    Mdouble q13 = 2 * q1 * q3;
    Mdouble q22 = q2 * q2;
    Mdouble q23 = 2 * q2 * q3;
    Mdouble q33 = q3 * q3;
    return Matrix3D(
            q00 + q11 - q22 - q33, q12 + q03, q13 - q02,
            q12 - q03, q00 - q11 + q22 - q33, q23 + q01,
            q13 + q02, q23 - q01, q00 - q11 - q22 + q33) * position;
}


/**
     * Calculates the distance from a wall through p0 whose normal is the vector (1,0,0).
     * Used for the calculation of the distance in InfiniteWalls
     */
Mdouble Quaternion::getDistance(const Vec3D p, const Vec3D p0) const
{
    const Mdouble q00 = q0 * q0;
    const Mdouble q02 = 2 * q0 * q2;
    const Mdouble q03 = 2 * q0 * q3;
    const Mdouble q11 = q1 * q1;
    const Mdouble q12 = 2 * q1 * q2;
    const Mdouble q13 = 2 * q1 * q3;
    const Mdouble q22 = q2 * q2;
    const Mdouble q33 = q3 * q3;
    return (q00 + q11 - q22 - q33) * (p0.X - p.X) +
           ((q12 + q03) * (p0.Y - p.Y) + (q13 - q02) * (p0.Z - p.Z));
}

/**
 * \todo move link to where it belongs
 * http://stackoverflow.com/questions/1171849/finding-quaternion-representing-the-rotation-from-one-vector-to-another
 * This is the same as rotateTensor, but than for MatrixSymmetric3D instead of SmallMatrix<3,3>.
 */
MatrixSymmetric3D Quaternion::rotateInverseInertiaTensor(const MatrixSymmetric3D& invI) const
{
    Mdouble a = 1 - 2 * q2 * q2 - 2 * q3 * q3;
    Mdouble b = 2 * q1 * q2 - 2 * q0 * q3;
    Mdouble c = 2 * q1 * q3 + 2 * q0 * q2;
    Mdouble d = 2 * q1 * q2 + 2 * q0 * q3;
    Mdouble e = 1 - 2 * q1 * q1 - 2 * q3 * q3;
    Mdouble f = 2 * q2 * q3 - 2 * q0 * q1;
    Mdouble g = 2 * q1 * q3 - 2 * q0 * q2;
    Mdouble h = 2 * q2 * q3 + 2 * q0 * q1;
    Mdouble i = 1 - 2 * q1 * q1 - 2 * q2 * q2;
    MatrixSymmetric3D ans = MatrixSymmetric3D(
            invI.XX * a * a + 2 * invI.XY * a * b + 2 * invI.XZ * a * c + invI.YY * b * b + 2 * invI.YZ * b * c +
            invI.ZZ * c * c,
            d * (invI.XX * a + invI.XY * b + invI.XZ * c) + e * (invI.XY * a + invI.YY * b + invI.YZ * c) +
            f * (invI.XZ * a + invI.YZ * b + invI.ZZ * c),
            g * (invI.XX * a + invI.XY * b + invI.XZ * c) + h * (invI.XY * a + invI.YY * b + invI.YZ * c) +
            i * (invI.XZ * a + invI.YZ * b + invI.ZZ * c),
            invI.XX * d * d + 2 * invI.XY * d * e + 2 * invI.XZ * d * f + invI.YY * e * e + 2 * invI.YZ * e * f +
            invI.ZZ * f * f,
            g * (invI.XX * d + invI.XY * e + invI.XZ * f) + h * (invI.XY * d + invI.YY * e + invI.YZ * f) +
            i * (invI.XZ * d + invI.YZ * e + invI.ZZ * f),
            invI.XX * g * g + 2 * invI.XY * g * h + 2 * invI.XZ * g * i + invI.YY * h * h + 2 * invI.YZ * h * i +
            invI.ZZ * i * i
    );
    return ans;
}

void Quaternion::rotateTensor(SmallMatrix<3, 3> I) const
{
    SmallMatrix<3, 3> A;
    getRotationMatrix(A);
    I = A * I * A.transpose();
}
