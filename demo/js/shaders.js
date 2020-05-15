window.shaders = window.shaders || {};


window.shaders.advectionShader = `
uniform vec2 res; //The width and height of our screen
uniform sampler2D bufferTexture; // The texture to advect
uniform sampler2D velocity; // Texture of 2D velocity vectors
uniform float dt; // current change in time;
uniform float dx; // grid spacing for square cells (distance in pixels between grid cells?)

vec4 lerp(vec4 v0, vec4 v1, float x);

// Linear interpolation function
vec4 lerp(vec4 v0, vec4 v1, float x) {
    return (1. - x) * v0 + x * v1;
}

void main() {
    vec2 pixel = gl_FragCoord.xy / res.xy;

    // Calculate the position at the next timestep using the advection equation
    vec4 vel = dt * texture2D(velocity, pixel);
    vec2 newPos = gl_FragCoord.xy - vel.xy;

    // Bounding the new position to be within the bounds (0.5, 0.5) and (res.x - 0.5, res.y - 0.5) so that its grid coordinates are between (0, 0) and (res.x - 1, res.y - 1).
    if (newPos.x < 0.5) {
    	newPos.x = 0.5;
    } else if (newPos.x > res.x - 0.5) {
    	newPos.x = res.x - 0.5;
    }
    if (newPos.y < 0.5) {
    	newPos.y = 0.5;
    } else if (newPos.y > res.y - 0.5) {
    	newPos.y = res.y - 0.5;
    }

    // The pixel at (0, 0) on the grid is located at (0.5, 0.5) on the screen
    // When interpolating, we want it to be located at (0, 0) instead
    // That way, a location halfway between pixel (0, 0) and pixel (1, 0) is located at (0.5, 0) instead of at (1.0, 0)
    // So we can use loc.x - floor(loc.x) as the x-value for lerp
    vec2 adjustedPos = newPos - 0.5;
    float x = adjustedPos.x - floor(adjustedPos.x);
    float y = adjustedPos.y - floor(adjustedPos.y);

    // Need to readd the 0.5 when sampling from the texture, since the pixel at (0, 0) on the grid is located at (0.5, 0.5) on the screen
    vec4 botLeft = texture2D(bufferTexture, (floor(adjustedPos) + 0.5) / res.xy);
    vec4 botRight = texture2D(bufferTexture, vec2(ceil(adjustedPos.x) + 0.5, floor(adjustedPos.y) + 0.5) / res.xy);
    vec4 topLeft = texture2D(bufferTexture, vec2(floor(adjustedPos.x) + 0.5, ceil(adjustedPos.y) + 0.5) / res.xy);
    vec4 topRight = texture2D(bufferTexture, (ceil(adjustedPos) + 0.5) / res.xy);

    // We manually interpolate instead of using WebGL's built in interpolation, because the built in interpolation apparently uses 24.8 precision texture interpolators
    // This means that there are only 256 maximum possible values between two adjacent pixels of a texture
    // Manually interpolating thus gives us more precision, reducing aliasing
    vec4 bot = lerp(botLeft, botRight, x);
    vec4 top =  lerp(topLeft, topRight, x);

    gl_FragColor = vec4(lerp(bot, top, y).xyz, 1.0);
}
`;

window.shaders.buoyancyShader = `

uniform vec2 res; //The width and height of our screen
uniform sampler2D temperature; // Texture of scalar temperatures
uniform sampler2D velocity; // Texture of 2D velocity vectors
uniform sampler2D density; // Texture of density
uniform float dt; // current change in time
uniform float T0; // ambient temperature
uniform float sigma; // buoyancy coefficient
uniform float kappa; // smoke weight


//outputs the value of the buoyant force which is only in the +y direction
void main() {
	vec2 pixel = gl_FragCoord.xy / res.xy;

	// CHANGE THIS IF YOU NO LONGER STORE TEMPERATURE IN THE X/R COMPONENT OF THE VECTOR
	float temp = texture2D(temperature, pixel).x;
	vec2 vel = texture2D(velocity, pixel).xy;
	// CHANGE THIS IF YOU NO LONGER STORE DENSITY IN THE Z/B COMPONENT OF THE VECTOR
	float dens = texture2D(density, pixel).z;

	gl_FragColor.xy = vel;
	// Buoyancy equation
	gl_FragColor.y += dt * (-1. * kappa * dens + sigma * (temp - T0)); // updates the velocity at pixel
	gl_FragColor.w = 1.;
}
`;

window.shaders.jacobiShader = `
uniform vec2 res; //The width and height of our screen
uniform sampler2D x; // The texture we want to calculate
uniform sampler2D b; // The texture we are comparing it against
uniform float alpha; // Alpha in the Jacobi equation
uniform float rBeta; // reciprocal of Beta in the Jacobi equation

void main() {
	vec2 pixel = gl_FragCoord.xy / res.xy;

	float xPixel = 1.0/res.x;//The size of a single pixel
	float yPixel = 1.0/res.y;

	vec4 right = texture2D(x, vec2(pixel.x + xPixel, pixel.y));
	vec4 left = texture2D(x, vec2(pixel.x - xPixel, pixel.y));
	vec4 up = texture2D(x, vec2(pixel.x, pixel.y + yPixel));
	vec4 down = texture2D(x, vec2(pixel.x, pixel.y - yPixel));

	vec4 bC = texture2D(b, pixel);

	// Jacobi equation
	gl_FragColor = (left + right + down + up + alpha * bC) * rBeta;
	gl_FragColor.w = 1.;
}
`;

window.shaders.divergenceShader = `

uniform vec2 res; // width and height of our screen
uniform sampler2D w; // Texture of 2D velocity vectors
uniform float dx; // Grid spacing for square grid cells

void main() {
	vec2 pixel = gl_FragCoord.xy / res.xy;

	float xPixel = 1.0/res.x; //The size of a single pixel
	float yPixel = 1.0/res.y;

	vec4 right = texture2D(w, vec2(pixel.x + xPixel, pixel.y));
	vec4 left = texture2D(w, vec2(pixel.x - xPixel, pixel.y));
	vec4 up = texture2D(w, vec2(pixel.x, pixel.y + yPixel));
	vec4 down = texture2D(w, vec2(pixel.x, pixel.y - yPixel));

	// Divergence equation
	gl_FragColor.r = 0.5 / dx * ((right.x - left.x) + (up.y - down.y));
	gl_FragColor.w = 1.;
}
`;

window.shaders.gradientShader = `

uniform vec2 res; // width and height of our screen
uniform sampler2D p; // Texture of scalar pressures. The x component stores the pressure
uniform sampler2D v; // Texture of 2D velocity vectors
uniform float dx; // Grid spacing for square cells
void main() {
	vec2 pixel = gl_FragCoord.xy / res.xy;

	float xPixel = 1.0/res.x;//The size of a single pixel
	float yPixel = 1.0/res.y;

	vec4 right = texture2D(p, vec2(pixel.x + xPixel, pixel.y));
	vec4 left = texture2D(p, vec2(pixel.x - xPixel, pixel.y));
	vec4 up = texture2D(p, vec2(pixel.x, pixel.y + yPixel));
	vec4 down = texture2D(p, vec2(pixel.x, pixel.y - yPixel));

	vec4 uNew = texture2D(v, pixel);
	// Gradient equation
	uNew.xy -= 0.5 / dx * vec2(right.x - left.x, up.x - down.x);

	gl_FragColor = uNew;
	gl_FragColor.w = 1.;
}
`;

window.shaders.createBoundaryShader = `
uniform vec2 res; // width and height of our screen
uniform sampler2D boundary; // Texture of floats. 1 means it's a boundary cell, 0 means it's not
uniform vec2 start; // The starting x and y coordinates of the bounding box of the new boundary
uniform vec2 end; // The ending x and y coordinates of the bounding box of the new boundary

void main() {
	vec2 pixel = gl_FragCoord.xy / res.xy;
	vec2 boundStart = start / res.xy;
	vec2 boundEnd = end / res.xy;

	// If the pixel is within the bounding box, create a new boundary cell
	if (pixel.y >= boundStart.y && pixel.y <= boundEnd.y && pixel.x >= boundStart.x && pixel.x <= boundEnd.x) {
		gl_FragColor.r = 1.;
	} else {
		// Otherwise use existing value
		gl_FragColor.r = texture2D(boundary, pixel).r;
	}
	gl_FragColor.w = 1.;
}
`

window.shaders.boundaryShader = `
uniform vec2 res; // width and height of our screen
uniform sampler2D texture; // The texture we are enforcing boundary conditions on
uniform sampler2D boundary; // Texture of floats. 1 means it's a boundary cell, 0 means it's not
uniform float mode; // 0 if continuity, 1 if opposite
uniform float scale; // -1 if texture if velocity, 1 if texture is pressure, temperature or density
void main() {
	vec2 pixel = gl_FragCoord.xy / res.xy;

	float xPixel = 1.0/res.x;//The size of a single pixel
	float yPixel = 1.0/res.y;

	vec2 left = vec2(pixel.x - xPixel, pixel.y);
	vec2 right = vec2(pixel.x + xPixel, pixel.y);
	vec2 top = vec2(pixel.x, pixel.y + yPixel);
	vec2 bot = vec2(pixel.x, pixel.y - yPixel);

	// If this cell is a boundary cell
	if (texture2D(boundary, pixel).x > 0.1) {
		// If the cell to the right of this cell is not a boundary cell
		if (texture2D(boundary, right).x < 0.1) {
			// If we want this cell to have the opposite horizontal velocity as its neighbor
			if (mode > 0.5) {
				gl_FragColor.x = scale * texture2D(texture, right).x;
				gl_FragColor.yz = texture2D(texture, right).yz;
			} else {
				// Otherwise, ensure continuity
				gl_FragColor.xyz = scale * texture2D(texture, right).xyz;
			}
		} else if (texture2D(boundary, left).x < 0.1) { // If the cell to the left of this cell is not a boundary cell
			// If we want this cell to have the opposite horizontal velocity as its neighbor
			if (mode > 0.5) {
				gl_FragColor.x = scale * texture2D(texture, left).x;
				gl_FragColor.yz = texture2D(texture, left).yz;
			} else {
				// Otherwise, ensure continuity
				gl_FragColor.xyz = scale * texture2D(texture, left).xyz;
			}
		} else if (texture2D(boundary, top).x < 0.1) { // If the cell above this cell is not a boundary cell
			// If we want this cell to have the opposite vertical velocity as its neighbor
			if (mode > 0.5) {
				gl_FragColor.y = scale * texture2D(texture, top).y;
				gl_FragColor.xz = texture2D(texture, top).xz;
			} else {
				// Otherwise, ensure continuity
				gl_FragColor.xyz = scale * texture2D(texture, top).xyz;
			}
		} else if (texture2D(boundary, bot).x < 0.1) { // If the cell below this cell is not a boundary cell
			// If we want this cell to have the opposite vertical velocity as its neighbor
			if (mode > 0.5) {
				gl_FragColor.y = scale * texture2D(texture, bot).y;
				gl_FragColor.xz = texture2D(texture, bot).xz;
			} else{
				// Otherwise, ensure continuity
				gl_FragColor.xyz = scale * texture2D(texture, bot).xyz;
			}
		} else { // If all of its neighbors are boundary cells, use the average of its neighbors
			gl_FragColor.xyz = scale * (texture2D(texture, right).xyz + texture2D(texture, left).xyz + texture2D(texture, top).xyz + texture2D(texture, bot).xyz) / 4.;
		}
	} else { // If this cell is not a boundary cell, use the current value
		gl_FragColor = texture2D(texture, pixel);
	}
	gl_FragColor.w = 1.;
}
`;

window.shaders.exForceDensityShader = `
uniform vec2 res; // width and height of our screen
uniform sampler2D bufferTexture; // Density texture
uniform sampler2D boundary; // Texture of floats. 1 means it's a boundary cell, 0 means it's not.
uniform vec3 smokeSource; // The x and y components contain the location of the smoke source. The z component contains the density.
uniform float radius; // desired impulse radius

void main() {
	vec2 pixel = gl_FragCoord.xy / res.xy;
	gl_FragColor = texture2D(bufferTexture, pixel);

	// Only add density if it's not a boundary cell
	if (texture2D(boundary, pixel).x < 0.1) {
		//Get the distance of the current pixel from the smoke source
		float dist = distance(smokeSource.xy,gl_FragCoord.xy);

		// Add to current density and make smoke blue.
		gl_FragColor.rgb += smokeSource.z * max(radius-dist,0.0);
		gl_FragColor.w = 1.;
	}
}
`;

window.shaders.exForceTemperatureShader = `
uniform vec2 res; // width and height of our screen
uniform sampler2D bufferTexture; // Temperature texture
uniform sampler2D boundary; // Texture of floats. 1 means it's a boundary cell, 0 means it's not.
uniform vec3 smokeSource; // The x and y components contain the location of the smoke source. The z component contains the temperature.
uniform float radius; // desired impulse radius

void main() {
	vec2 pixel = gl_FragCoord.xy / res.xy;
	gl_FragColor = texture2D(bufferTexture, pixel);

	// Only add temperature if it's not a boundary cell
	if (texture2D(boundary, pixel).x < 0.1) {
		//Get the distance of the current pixel from the smoke source
		float dist = distance(smokeSource.xy,gl_FragCoord.xy);

		// Add to current temperature and make smoke red
		gl_FragColor.rgb += smokeSource.z * max(radius-dist,0.0);
		gl_FragColor.w = 1.;
	}
}
`;



window.shaders.exForceVelocityShader = `
uniform vec2 res; // width and height of our screen
uniform sampler2D bufferTexture; // Velocity texture
uniform sampler2D boundary; // Texture of floats. 1 means it's a boundary cell, 0 means it's not.
uniform vec3 smokeSource; // The x and y components contain the location of the smoke source. The z component contains the density.
uniform float radius; // desired impulse radius

void main() {
	vec2 pixel = gl_FragCoord.xy / res.xy;
	gl_FragColor = texture2D(bufferTexture, pixel);

	// Only add velocity if the cell is not a boundary cell
	if (texture2D(boundary, pixel).x < 0.1) {
		//Get the distance of the current pixel from the smoke source
		float dist = distance(smokeSource.xy,gl_FragCoord.xy);
		vec2 direction = normalize(gl_FragCoord.xy - smokeSource.xy);

		// Add to current velocity
		gl_FragColor.rg += smokeSource.z * max(radius-dist,0.0) * direction;
		gl_FragColor.w = 1.;
	}
}
`;

window.shaders.curlShader = `
uniform vec2 res; // width and height of our screen
uniform sampler2D u; // Texture of velocity
uniform float dx; // Grid spacing for square cells
void main() {
	vec2 pixel = gl_FragCoord.xy / res.xy;

	float xPixel = 1.0/res.x;//The size of a single pixel
	float yPixel = 1.0/res.y;

	gl_FragColor = vec4(0, 0, 0, 1);

	vec4 right = texture2D(u, vec2(pixel.x + xPixel, pixel.y));
	vec4 left = texture2D(u, vec2(pixel.x - xPixel, pixel.y));
	vec4 up = texture2D(u, vec2(pixel.x, pixel.y + yPixel));
	vec4 down = texture2D(u, vec2(pixel.x, pixel.y - yPixel));

	// 2D cross product equation
	gl_FragColor.z = 0.5 / dx * (right.y - left.y - (up.x - down.x));
}
`;

window.shaders.normalizeShader = `
uniform vec2 res; // width and height of our screen
uniform sampler2D w; // Texture of vorticity
uniform float dx; // Grid spacing for square cells
void main() {
	vec2 pixel = gl_FragCoord.xy / res.xy;

	float xPixel = 1.0/res.x;//The size of a single pixel
	float yPixel = 1.0/res.y;

	gl_FragColor = vec4(0, 0, 0, 1);

	vec4 right = texture2D(w, vec2(pixel.x + xPixel, pixel.y));
	vec4 left = texture2D(w, vec2(pixel.x - xPixel, pixel.y));
	vec4 up = texture2D(w, vec2(pixel.x, pixel.y + yPixel));
	vec4 down = texture2D(w, vec2(pixel.x, pixel.y - yPixel));

	// Gradient equation
	vec2 temp = 0.5 / dx * vec2(abs(right.z) - abs(left.z), abs(up.z) - abs(down.z));
	// Avoid division by 0
	if (temp.x > 0.000001) {
		gl_FragColor.r = temp.x / abs(temp.x);
	}
	if (temp.y > 0.000001) {
		gl_FragColor.g = temp.y / abs(temp.y);
	}
}
`;

window.shaders.vorticityForceShader = `
uniform vec2 res; // width and height of our screen
uniform sampler2D w; // Texture of vorticity
uniform sampler2D psi; // Texture of normalized vorticity
uniform sampler2D velocity; // Texture of velocity;
uniform float dx; // Grid spacing for square cells
uniform float epsilon; // Scale parameter
uniform float dt; // Timestep size
void main() {
	vec2 pixel = gl_FragCoord.xy / res.xy;

	float xPixel = 1.0/res.x;//The size of a single pixel
	float yPixel = 1.0/res.y;

	gl_FragColor = texture2D(velocity, pixel);

	vec4 a = texture2D(psi, pixel);
	vec4 b = texture2D(w, pixel);

	// Cross product equation
	gl_FragColor.x += epsilon * dx * dt * (a.y * b.z - a.z * b.y);
	gl_FragColor.y += epsilon * dx * dt * (a.z * b.x - a.x * b.z);
	gl_FragColor.z += epsilon * dx * dt * (a.x * b.y - a.y * b.x);
}
`;

window.shaders.displayShader = `
uniform vec2 res; // width and height of our screen
uniform sampler2D bufferTexture; // Texture to draw
uniform sampler2D boundary; // Boundary texture so we can display boundaries as a separate color
uniform float displayBoundaries; // 1 if we should show boundaries in gray, 0 otherwise.
uniform vec3 color; // Rgb scale parameters
uniform float scaleColor; // Scale factor to make sure we don't saturate, but also make sure we can see the texture values
uniform float biasColor; // Bias term to make sure we can render negative texture values
uniform float showInsideBoundaries; // 1 if we should show the values of buffertexture inside the boundary, 0 otherwise
void main() {
	vec2 pixel = gl_FragCoord.xy / res.xy;

	// If this cell is not a boundary cell or if we shouldn't show boundaries, then use the existing value as the color.
	if (texture2D(boundary, pixel).x < 0.1) {
		gl_FragColor.rgb = texture2D(bufferTexture, pixel).rgb * color * scaleColor + biasColor;
		gl_FragColor.w = 1.;
	} else if (displayBoundaries > 0.1) {
		// Show boundaries in gray.
		gl_FragColor = vec4(0.25, 0.25, 0.25, 1);
	} else {
		if (showInsideBoundaries > 0.1) {
			gl_FragColor.rgb = texture2D(bufferTexture, pixel).rgb * color * scaleColor + biasColor;
		} else {
			gl_FragColor = vec4(biasColor, biasColor, biasColor, 1);
		}
	}
}
`

window.shaders.resetShader = `
uniform vec2 res; // width and height of our screen
uniform sampler2D bufferTexture; // Texture to remove boundary conditions from
uniform sampler2D boundaryA; // Previous boundary texture
uniform sampler2D boundaryB; // Boundary texture after reset
void main() {
	vec2 pixel = gl_FragCoord.xy / res.xy;

	// If the boundary was removed, reset the value in the texture
	if (texture2D(boundaryA, pixel).x > 0.1 && texture2D(boundaryB, pixel).x < 0.1) {
		gl_FragColor = vec4(0, 0, 0, 1);
	} else {
		// Otherwise, use the current value
		gl_FragColor = texture2D(bufferTexture, pixel);
	}
}
`
