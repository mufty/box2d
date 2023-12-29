/*
 * Copyright (c) 2006-2007 Erin Catto http://www.gphysics.com
 *
 * This software is provided 'as-is', without any express or implied
 * warranty.  In no event will the authors be held liable for any damages
 * arising from the use of this software.
 * Permission is granted to anyone to use this software for any purpose,
 * including commercial applications, and to alter it and redistribute it
 * freely, subject to the following restrictions:
 * 1. The origin of this software must not be misrepresented; you must not
 * claim that you wrote the original software. If you use this software
 * in a product, an acknowledgment in the product documentation would be
 * appreciated but is not required.
 * 2. Altered source versions must be plainly marked as such, and must not be
 * misrepresented as being the original software.
 * 3. This notice may not be removed or altered from any source distribution.
 */

package box2D.dynamics.contacts;

import box2D.collision.B2ManifoldType;
import box2D.collision.B2Collision;
import box2D.common.math.B2Vec2;
import box2D.common.math.B2Math;
import box2D.collision.B2Manifold;
import box2D.collision.shapes.B2CircleShape;
import box2D.collision.shapes.B2EdgeShape;
import box2D.collision.shapes.B2ShapeType;
import box2D.common.B2Settings;
import box2D.common.math.B2Transform;
import box2D.dynamics.B2Body;
import box2D.dynamics.B2Fixture;
import box2D.dynamics.contacts.B2Contact;

/**
 * @private
 */
class B2EdgeAndCircleContact extends B2Contact
{
	static public function create(allocator:Dynamic):B2Contact
	{
		return new B2EdgeAndCircleContact();
	}

	static public function destroy(contact:B2Contact, allocator:Dynamic):Void
	{
		//
	}

	public override function reset(fixtureA:B2Fixture = null, fixtureB:B2Fixture = null):Void
	{
		super.reset(fixtureA, fixtureB);
		B2Settings.b2Assert(fixtureA == null || fixtureA.getType() == B2ShapeType.EDGE_SHAPE);
		B2Settings.b2Assert(fixtureB == null || fixtureB.getType() == B2ShapeType.CIRCLE_SHAPE);
	}

	// ~b2EdgeAndCircleContact() {}

	public override function evaluate():Void
	{
		var bA:B2Body = m_fixtureA.getBody();
		var bB:B2Body = m_fixtureB.getBody();
		b2CollideEdgeAndCircle(m_manifold, cast(m_fixtureA.getShape(), B2EdgeShape), bA.m_xf, cast(m_fixtureB.getShape(), B2CircleShape), bB.m_xf);
	}

	private function b2CollideEdgeAndCircle(manifold:B2Manifold, edgeA:B2EdgeShape, xfA:B2Transform, circleB:B2CircleShape, xfB:B2Transform):Void
	{
        manifold.m_pointCount = 0;

        // Compute circle in frame of edge
        var Q:B2Vec2 = B2Math.mulX(xfA, B2Math.mulX(xfB, circleB.m_p));

        var A:B2Vec2 = edgeA.getVertex1(), B:B2Vec2 = edgeA.getVertex2();
        var e:B2Vec2 = B2Math.subtractVV(B, A); //b2Vec2 e = B - A;

        // Normal points to the right for a CCW winding
        var n:B2Vec2 = new B2Vec2(e.y, -e.x);
        var offset:Float = B2Math.dot(n, B2Math.subtractVV(Q, A)); //float offset = b2Dot(n, Q - A);

        var oneSided = false; //TODO this isn't supported in B2EdgeShape // bool oneSided = edgeA->m_oneSided;
        if (oneSided && offset < 0.0)
        {
            return;
        }

        // Barycentric coordinates
        var u:Float = B2Math.dot(e, B2Math.subtractVV(B, Q)); //float u = b2Dot(e, B - Q);
		var v:Float = B2Math.dot(e, B2Math.subtractVV(Q, A)); //float v = b2Dot(e, Q - A);

        var radius:Float = edgeA.m_radius + circleB.m_radius;

        //NOTE we don't need this?
        /*var cf:B2ContactFeature = {
            indexB: 0,
            typeB: B2ContactFeatureType.e_vertex,
            indexA: 0,
            typeA: B2ContactFeatureType.e_vertex
        };*/
		/*
		b2ContactFeature cf;
		cf.indexB = 0;
		cf.typeB = b2ContactFeature::e_vertex;
        */
		
        // Region A
        if (v <= 0.0)
        {
            var P:B2Vec2 = A;
            var d:B2Vec2 = B2Math.subtractVV(Q, P);
            var dd:Float = B2Math.dot(d, d);
            if(dd > radius * radius){
                return;
            }

            // Is there an edge connected to A?
            //TODO onesided not supported
            /*
            if (edgeA->m_oneSided)
			{
				b2Vec2 A1 = edgeA->m_vertex0;
				b2Vec2 B1 = A;
				b2Vec2 e1 = B1 - A1;
				float u1 = b2Dot(e1, B1 - Q);
				
				// Is the circle in Region AB of the previous edge?
				if (u1 > 0.0f)
				{
					return;
				}
			}
            */

            /*cf.indexA = 0;
            cf.typeA = B2ContactFeatureType.e_vertex;*/
            manifold.m_pointCount = 1;
            manifold.m_type = B2ManifoldType.CIRCLES;
            manifold.m_localPlaneNormal.setZero(); //Note not sure if fine it could be null
            manifold.m_localPoint = P;
            manifold.m_points[0].m_id._key = 0; //Note not exactly null safe
            //manifold.m_points[0].m_id.cf = cf;
            manifold.m_points[0].m_localPoint = circleB.m_p;
            return;
        }

        // Region B
        if(u <= 0.0){
            var P:B2Vec2 = B;
            var d:B2Vec2 = B2Math.subtractVV(Q, P);
            var dd:Float = B2Math.dot(d, d);

            if (dd > radius * radius) {
                return;
            }

            // Is there an edge connected to B?
            //TODO onesided not supported
            /*if (edgeA->m_oneSided)
            {
                b2Vec2 B2 = edgeA->m_vertex3;
                b2Vec2 A2 = B;
                b2Vec2 e2 = B2 - A2;
                float v2 = b2Dot(e2, Q - A2);
                
                // Is the circle in Region AB of the next edge?
                if (v2 > 0.0f)
                {
                    return;
                }
            }*/

            /*
            cf.indexA = 1;
			cf.typeA = b2ContactFeature::e_vertex;
            */
            manifold.m_pointCount = 1;
			manifold.m_type = B2ManifoldType.CIRCLES;
			manifold.m_localPlaneNormal.setZero(); //Note not sure if fine it could be null
			manifold.m_localPoint = P;
			manifold.m_points[0].m_id._key = 0; //Note not exactly null safe
            //manifold.m_points[0].m_id.cf = cf;
            manifold.m_points[0].m_localPoint = circleB.m_p;
			return;
        }

        // Region AB
        var den:Float = B2Math.dot(e, e);
        B2Settings.b2Assert(den > 0.0);
        var uA = B2Math.mulFV(u, A);
        var vB = B2Math.mulFV(v, B);
        var uAvB = B2Math.addVV(uA, vB);
        var P:B2Vec2 = B2Math.mulFV((1.0 / den), uAvB); //b2Vec2 P = (1.0f / den) * (u * A + v * B);
        var d:B2Vec2 = B2Math.subtractVV(Q, P);
        var dd:Float = B2Math.dot(d, d);
        if (dd > radius * radius)
        {
            return;
        }

        if (offset < 0.0)
        {
            n.set(-n.x, -n.y);
        }
        n.normalize();
        /*	
		cf.indexA = 0;
		cf.typeA = b2ContactFeature::e_face;
        */
		manifold.m_pointCount = 1;
		manifold.m_type = B2ManifoldType.FACE_A;
		manifold.m_localPlaneNormal = n;
		manifold.m_localPoint = A;
		manifold.m_points[0].m_id._key = 0;
		//manifold->points[0].id.cf = cf;
		manifold.m_points[0].m_localPoint = circleB.m_p;
	}
}
