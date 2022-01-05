#include <limits>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <SDL.h>
#include <SDL_image.h>
#undef main

#include <cmath>
const int SCREEN_WIDTH = 1024;
const int SCREEN_HEIGHT = 768;
const int centerX = SCREEN_WIDTH / 2;
const int centerY = SCREEN_HEIGHT / 2;
const double pi = 3.1415926535;



SDL_Window* Window = NULL;
SDL_Surface* WindowSurface = NULL;
SDL_Surface* WorkingSurface = NULL;

unsigned char rawBuffer[4 * SCREEN_WIDTH * SCREEN_HEIGHT];


//ANY SORT OF DOT PRODUCT NEEDS AN ORIGIN TO BE ADDED TO IT BRUH DUH


struct Vec3f {
    public:
        float x, y, z;
        float registers[3];
        Vec3f(const float& i = 0, const float& j = 0, const float& k = 0) {
            x = i;
            y = j;
            z = k;

            registers[0] = x;
            registers[1] = y;
            registers[2] = z;
        }
        float magnitude() {
            return pow(x * x + y * y + z * z, .5);
        }
        //outputs string representation of vector, with newline character after.
        std::string toString() {
            return "<" + std::to_string(x) + "," + std::to_string(y) + "," + std::to_string(z) + ">\n";
        }
        float operator [](int i) const { return registers[i]; }
        float& operator [](int i) { return registers[i]; }
        friend Vec3f operator+(const Vec3f& a, const Vec3f& b) {
            return Vec3f(a.x + b.x, a.y + b.y, a.z + b.z);
        }
        friend Vec3f operator-(const Vec3f& a, const Vec3f& b) {
            return Vec3f(a.x - b.x, a.y - b.y, a.z - b.z);
        }
        
        friend Vec3f operator*(const double& a, const Vec3f& b) {
            return Vec3f(a*b.x, a*b.y, a*b.z);
        }
        friend Vec3f operator*(const Vec3f& b, const double& a) {
            return Vec3f(a * b.x, a * b.y, a * b.z);
        }
        
        friend float operator*(const Vec3f& a, const Vec3f& b) {
            return a.x * b.x + a.y * b.y + a.z * b.z;
        }
};

//temp object to see if idea works
Vec3f rayProjections[SCREEN_WIDTH * SCREEN_HEIGHT];
//temp light source to see if idea works
Vec3f tempLightSource = Vec3f(1/sqrt(3),1/sqrt(3),-1/sqrt(3));
//global origin to save my ass this my last resort (havent used it yet hmm)
Vec3f globalOrigin = Vec3f(0, 0, 0);


//I DID ITWOOOO

struct Sphere {
public:
    Vec3f center;
    float radius;
    Sphere(const float& i = 0, const float& j = 0, const float& k = 0, const float& rad = 0) {
        center = Vec3f(i, j, k);
        radius = rad;
    }
    //direction is a unit vector always
    //distance is the distance to the first collision IF TRUE
    //but im going to use distance to create shadows rn
    bool ray_intersection(Vec3f& origin, Vec3f& direction, float& dist) {

        //REPEATED CALCULATIONS HERE CAN BE OPTIMIZED


        Vec3f B_A = center - origin;
        if ((B_A * direction) < 0) {
            return false;
        }
        Vec3f projection = origin + (B_A * direction)*direction;
        Vec3f dis = projection - center;


        dist = dis.magnitude();
        if (dist <= radius){
            

            //I NEEDED TO SUBTRACT THE ORIGIN FROM PROJECTION TO KEEP PROJECTION IN THE RELATIVE OBJECT-CAMERA SPACE
            float scaleFactor = (projection - origin).magnitude();
            float subtract = radius * radius - dis.magnitude() * dis.magnitude();
            subtract = pow(subtract, .5);
            //should be magnitude of vector pointing to first contact now
            scaleFactor -= subtract;
            //now vector pointing to first contact INEEEDED AN ORGINI ADDED HERE OH MY LORD IT MAKES SO MUCH SENSE NOW
            //nope better but still wrong, but i do know first contact should proportionally decrease as I get closer
            //FIXED IT BY SUBTRACTING ORIGINWOOO
            Vec3f firstContact =  origin + scaleFactor * direction;
            Vec3f normalVector = firstContact - center;
            normalVector = normalVector * (1.0 / normalVector.magnitude());
            dist =  normalVector * tempLightSource;
         

            if (dist < 0) {
                dist = abs(dist);
            }
            else {
                dist = 0;
            }

            
            //if i could just find whatever transformations needed to maintain the square on the surface of the sphere than i have won
            //I HAVE WON
            if (abs(firstContact[0]) < .1 && abs(firstContact[1] - 25) < .1 && abs(firstContact[2]) < .1) {
                //std::cout << (firstContact + origin).magnitude() << "," << firstContact.magnitude() << " " << projection.toString();
                dist = 0;
            }

            return true;
        }
        return false;
    }
};




uint32_t pack_color(const uint8_t r, const uint8_t g, const uint8_t b, uint8_t a) {
    return (a << 24) + (b << 16) + (g << 8) + r;
}
void unpack_color(const uint32_t color, uint8_t& r, uint8_t& g, uint8_t& b, uint8_t& a) {
    r = (color >> 0);
    g = (color >> 8);
    b = (color >> 16);
    a = (color >> 24);
}

void renderImage() {
    std::vector<Vec3f> framebuffer(SCREEN_WIDTH * SCREEN_HEIGHT);

    for (size_t j = 0; j < SCREEN_HEIGHT; j++) {
        for (size_t i = 0; i < SCREEN_WIDTH; i++) {
            framebuffer[i + j * SCREEN_WIDTH] = Vec3f(j / float(SCREEN_HEIGHT), i / float(SCREEN_WIDTH), 0);
        }
    }
    std::ofstream ofs; // save the framebuffer to file
    ofs.open("./out.ppm", std::ofstream::out | std::ofstream::binary);
    ofs << "P6\n" << SCREEN_WIDTH << " " << SCREEN_HEIGHT << "\n255\n";
    for (size_t i = 0; i < SCREEN_HEIGHT * SCREEN_WIDTH; ++i) {
        for (size_t j = 0; j < 3; j++) {
            ofs << (char)(255 * std::max(0.f, std::min(1.f, framebuffer[i][j])));
        }
    }
    ofs.close();
}
void fillRawBuffer(unsigned char* rawBuffer, std::vector<Vec3f>& frameBuffer) {
    float scale = 255;
    for (size_t i = 0; i < frameBuffer.size(); i++) {
        rawBuffer[i * 4] = (uint8_t)(scale * frameBuffer[i][0]);
        rawBuffer[i * 4 + 1] = (uint8_t)(scale * frameBuffer[i][1]);
        rawBuffer[i * 4 + 2] = (uint8_t)(scale * frameBuffer[i][2]);
        //for now a is always 0
        rawBuffer[i * 4 + 3] = 0;
    }
}

void fillRawBuffer(unsigned char* rawBuffer) {
    for (size_t i = 0; i < SCREEN_HEIGHT*SCREEN_WIDTH; i++) {
        rawBuffer[i * 4] = 255;
        rawBuffer[i * 4 + 1] = 255;
        rawBuffer[i * 4 + 2] = 255;
        //for now a is always 0
        rawBuffer[i * 4 + 3] = 0;
    }
}
void renderFrame(std::vector<Vec3f> &framebuffer) {
    for (size_t j = 0; j < SCREEN_HEIGHT; j++) {
        for (size_t i = 0; i < SCREEN_WIDTH; i++) {
            //ITS RGB INWARDSSSS
            //framebuffer[i + j * SCREEN_WIDTH] = Vec3f(0, 0, 1);
            float tempx = i / float(SCREEN_WIDTH);
            float tempy = j / float(SCREEN_HEIGHT);
            framebuffer[i + j * SCREEN_WIDTH] = Vec3f(0, tempx, tempy);
        }
    }
    fillRawBuffer(rawBuffer, framebuffer);
}

//camo is the camera position with respect to the screen to create the projection, camdir is its direction,
//dist is the distance from the camera to the screen
//function takes parameters into account and fills rayProjections with appropriate rays
void castRays(Vec3f camdir, float dist) {



    float theta;
    //theta and phi value conversions of the direction vector (need to make sure these conversions work for all theta n phi)
    if (camdir[0] == 0) {
        theta = pi/2.0;
    }
    else {
        theta = atan(camdir[1]/camdir[0]);
    }
    //using fact that rho is 1 (unit vector)
    float phi = acos(camdir[2]);

   

    //octant corrections for theta and phi below
    //also fixing for negative zero by adding + 0.0


    if (camdir[0] >= -0.0000001 && camdir[1] >= -0.0000001) {
        theta = abs(theta);
    }


    if (camdir[0] < -0.0000001 && camdir[1] >= 0) {
        std::cout << "first trigger \n";
        theta = pi - abs(theta);
    }

    else if (camdir[0] < -0.0000001 && camdir[1] < -0.0000001) {
        std::cout << "second trigger \n";
        theta = pi + theta;
    }
    else if (camdir[0] > 0 && camdir[1] < -0.0000001) {
        std::cout << "third trigger \n";
        theta = 2 * pi + theta;
    }

    if (camdir[2] <= 0) {
        phi = abs(phi);
    }







    Vec3f x_hat = Vec3f(-sin(theta), cos(theta), 0);
    Vec3f y_hat = Vec3f(-cos(phi) * cos(theta), -cos(phi) * sin(theta), sin(phi));

    Vec3f scaledX, scaledY,newDir;

    std::cout << "theta: " << theta << "\n";
    std::cout << "phi: " << phi << "\n";
    std::cout << "x_hat: " << x_hat.toString();
    std::cout << "y_hat: " << y_hat.toString();
    for (size_t i = 0; i < SCREEN_HEIGHT; i++) {
        for (size_t j = 0; j < SCREEN_WIDTH; j++) {
            //actually i think the diff values will flip signs, which will work in my favor in adding or subtracting 
                //we are on left half of projection
                    //we are on top half of projection
            
            //first get scale values 
            float scaleX = float(centerX - (int)j);
            float scaleY = float(centerY - (int)i);

            //then use defined unit vectors...
            scaledX = scaleX * x_hat;
            scaledY = scaleY * y_hat;
            
            //finally add to original vector to get new vector pointing in that direction
            newDir = camdir*dist + scaledX + scaledY;
            //scale the new direction to unit veector
            newDir = newDir * (1 / newDir.magnitude());
           
            rayProjections[j + i * SCREEN_WIDTH] = newDir;
        }
    }
}


void drawSphere(Sphere& circ, Vec3f& camo, float& dist, int* color, bool debug = false) {

if (debug){
    for (size_t j = 0; j < SCREEN_HEIGHT; j++) {
        for (size_t i = 0; i < SCREEN_WIDTH; i++) {
            if (circ.ray_intersection(camo, rayProjections[i + SCREEN_WIDTH * j], dist)) {

                if (dist <= 1) {
                    ;//do nothing
                }if (dist == 2) {
                    color[0] = color[1] = 0;
                    color[2] = 255;
                    dist = 1;
                }
                else if(dist == 3) {
                    color[0] = color[2] = 0;
                    color[1] = 255;
                    dist = 1;
                }
                else if(dist == 4) {
                    color[1] = color[2] = 0;
                    color[0] = 255;
                    dist = 1;
                }


                rawBuffer[(i + SCREEN_WIDTH * j) * 4] = color[0] * dist;
                rawBuffer[(i + SCREEN_WIDTH * j) * 4 + 1] = color[1] * dist;
                rawBuffer[(i + SCREEN_WIDTH * j) * 4 + 2] = color[2] * dist;
                rawBuffer[(i + SCREEN_WIDTH * j) * 4 + 3] = 255;

            }
        }
    }

}
else {
    for (size_t j = 0; j < SCREEN_HEIGHT; j++) {
        for (size_t i = 0; i < SCREEN_WIDTH; i++) {
            if (circ.ray_intersection(camo, rayProjections[i + SCREEN_WIDTH * j], dist)) {
                rawBuffer[(i + SCREEN_WIDTH * j) * 4] = color[0] * dist;
                rawBuffer[(i + SCREEN_WIDTH * j) * 4 + 1] = color[1] * dist;
                rawBuffer[(i + SCREEN_WIDTH * j) * 4 + 2] = color[2] * dist;
                rawBuffer[(i + SCREEN_WIDTH * j) * 4 + 3] = 255;

            }
        }
    }
}
}

int init() {
    //return flag
    bool success = true;

    //initialize the SDL
    if (SDL_Init(SDL_INIT_VIDEO) < 0) {
        printf("SDL could not initialize! SDL_Error: %s\n", SDL_GetError());
        success = false;
    }
    else {
        //successful init
        //create window
        Window = SDL_CreateWindow("Ray Tracer", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, SCREEN_WIDTH, SCREEN_HEIGHT, SDL_WINDOW_SHOWN);
        if (Window == NULL)
        {
            printf("Window could not be created! SDL_Error: %s\n", SDL_GetError());
            success = false;
        }
        else
        {
            //Get window surface
            WindowSurface = SDL_GetWindowSurface(Window);
        }
    }
    return success;
}


void close() {
    //Deallocate surface
    SDL_FreeSurface(WorkingSurface);
    WorkingSurface = NULL;

    //Destroy window
    SDL_DestroyWindow(Window);
    Window = NULL;

    //Quit SDL subsystems
    SDL_Quit();
}

//this function keeps a static global list of events for me to iterate thru
std::vector<SDL_Event>& GetFrameEvents()
{
    static std::vector<SDL_Event> frame_events;
    return frame_events;
}
//event handler
SDL_Event e;


int main() {

    /*
    
    
    
    MUST LIMIT PHI BETWEEN 0 AND PI AND THETA CAN RANGE FROM 0 TO 2PI
        
    
    */
    
    //circ.ray_intersection(origin, dir, dist);


    /*
    for (int i = 0; i < SCREEN_HEIGHT; i++) {
    std::cout << "i: " << i << "vec: " << rayProjections[SCREEN_WIDTH * (i)].toString();
    }
    */
    //starting sdl...

    if (!init()) {
        std::cout << "\n SDL INITIALIZE FAILED \n";
    }
    













    bool quit = false;

    float theta = pi/2;
    float phi = pi/2;
    float rho = 1;
    float x, y, z;
    x = y = z = 0;
    Vec3f dir;
    //main loop
    while (!quit) {

        //clearing buffer
        fillRawBuffer(rawBuffer);

        //looping through event queue
        while (SDL_PollEvent(&e) != 0) {
            //user requests quit so i dont even add to queue
            if (e.type == SDL_QUIT) {
                quit = true;
            }
            else {
                //where im adding to queue
                GetFrameEvents().push_back(e);
            }
        }



        /*
        
        TODO:
            GET MINECRAFT LIKE MOVEMENT THAT INCOPORATES THE MOUSE
            FIND OUT IF CAMERA ORIGINS ARE VALID FOR ALL OCTANTS (easy) AND FIX IF NOT (easyish)
            AMEND CAST RAYS TO BE ABLE TO USE MULTIPLE "CAMERA" SOURCES TO EVENTUALLY CREATE LIGHT SOURCES
                + CREATE LIGHT PROJECTIONS ON OBJECTS AS A RESULT
        
        */


        const Uint8* keyboard = SDL_GetKeyboardState(NULL);
        if (keyboard[SDL_SCANCODE_LEFT]) {
            theta += .1;
        }
        else if (keyboard[SDL_SCANCODE_RIGHT]) {
            theta -= .1;
        }
        if (keyboard[SDL_SCANCODE_UP]) {
            phi -= .1;
        }
        else if (keyboard[SDL_SCANCODE_DOWN]) {
            phi += .1;
        }
        //z is forward, x is backward
        if (keyboard[SDL_SCANCODE_Z]) {
            rho += 1;
        }
        else if (keyboard[SDL_SCANCODE_X]) {
            rho -= 1;
        }

        if (keyboard[SDL_SCANCODE_W]) {
            y += 1;
        }
        else if (keyboard[SDL_SCANCODE_S]) {
            y -= 1;
        }
        if (keyboard[SDL_SCANCODE_A]) {
            x -= 1;
        }
        else if (keyboard[SDL_SCANCODE_D]) {
            x += 1;
        }
        if (keyboard[SDL_SCANCODE_LSHIFT]) {
            z -= 1;
        }
        else if (keyboard[SDL_SCANCODE_SPACE]) {
            z += 1;
        }

        if (phi < 0) {
            phi = 0;
        }
        //avoid issues with unknown orientation
        else if (phi > pi - 0.00002) {
            phi = pi - 0.00002;
        }

        //std::vector<Vec3f>* frameBuffer = new std::vector<Vec3f>(SCREEN_WIDTH * SCREEN_HEIGHT);
//REMEBER RGB VALUES ARE REVERSED TO BGR
        int brown[3] = { 30, 105, 210 };
        int blue[3] = { 255, 0 ,0 };

        //renderFrame(*frameBuffer);
        Sphere circ = Sphere(0, 30, 0, 5);

        Sphere circ2 = Sphere(0, 100, 0, 20);

        //camera origin (seems to only matter in the intersection calculations)
        Vec3f origin = Vec3f(rho*sin(phi) * cos(theta) + x, rho*sin(phi) * sin(theta) + y, rho*cos(phi) + z);
        //camera pointing in user direction
        dir = Vec3f(sin(phi)*cos(theta),sin(phi)*sin(theta),cos(phi));
        std::cout << "DIRECTION: " << dir.toString();
        std::cout << "CAMORIGIN: " << origin.toString();
        float dist;

        //600 controls the distance between the screen and the camera, basically controlling for FOV
        castRays(dir, 600);

        //regular main loop calculations done every frame go here, such as drawing the frame

        

        //will have to learn how to stage objects soon
        drawSphere(circ2, origin, dist, brown);
        drawSphere(circ, origin, dist, blue, false);
   



        WorkingSurface = SDL_CreateRGBSurfaceFrom((void*)rawBuffer, SCREEN_WIDTH, SCREEN_HEIGHT, 32, 4 * SCREEN_WIDTH, 0, 0, 0, 0);
        SDL_BlitSurface(WorkingSurface, NULL, WindowSurface, NULL);
        SDL_UpdateWindowSurface(Window);
        //SDL_Delay(10000);
        //total events for the frame cleared here
        GetFrameEvents().clear();

        //capped to run at 60 frames per second
        SDL_Delay(16);
    }

    close();
    return 0;








}