

#include <time.h>
#include <stdlib.h>
#include <math.h>
#include "BandLimitedOsc.h"

#ifndef M_PI
#define M_PI (3.14159265358979323846264338327950288)
#endif

template<typename T>
inline int64_t bitwiseOrZero(const T &t) {
    return static_cast<int64_t>(t) | 0;
}

BandLimitedOsc::BandLimitedOsc() {
    m_wavetype = 2;
    m_freq = 1.f;
    m_sharp = 0.f;
    srand((unsigned int)time(0));
    m_sah_last_value = 0.f;
    m_sah_current_value = (rand() / (float)RAND_MAX) * 2.f - 1.f;
    phase1 = 0.f;
    phase2 = 0.f;
    pulse_width = 0.f;
    value1 = 0.f;
    value2 = 0.f;
}

void BandLimitedOsc::setup(float sampleRate) {
    m_sampleRate = sampleRate;
    m_sampleRate = 48000;
    m_oneOverSr = 1.f / m_sampleRate;
    m_twopi = 2.f * M_PI;
    m_oneOverPiOverTwo = 1.f / (M_PI / 2.f);
    m_oneOverPi = 1.f / M_PI;
    m_srOverFour = m_sampleRate / 4.f;
    m_srOverEight = m_sampleRate / 8.f;
    m_srOverTwo = m_sampleRate / 2.f;
    m_pointer_pos = m_sah_pointer_pos = 0.f;
}

BandLimitedOsc::~BandLimitedOsc() {}

void BandLimitedOsc::reset() {
    m_pointer_pos = m_sah_pointer_pos = 0.f;
}

float BandLimitedOsc::_clip(float x) {
    if (x < 0.f) {
        x += 1.f;
    } else if (x >= 1.f) {
        x -= 1.f;
    }
    return x;
}

void BandLimitedOsc::setWavetype(int type) {
    if (type != m_wavetype) {
        type = type < 0 ? 0 : type > 7 ? 7 : type;
        m_wavetype = type;
    }
}

void BandLimitedOsc::setFreq(float freq) {
    m_freq = freq < 0.00001f ? 0.00001f : freq > m_srOverFour ? m_srOverFour : freq;
}

void BandLimitedOsc::setSharp(float sharp) {
    m_sharp = sharp < 0.f ? 0.f : sharp > 1.f ? 1.f : sharp;
}

double BandLimitedOsc::poly_blep(double t)
{
    double dt = m_pointer_pos  / M_PI;
    // 0 <= t < 1
    if (t < dt) {
        t /= dt;
        return t+t - t*t - 1.0;
    }
    // -1 < t < 0
    else if (t > 1.0 - dt) {
        t = (t - 1.0) / dt;
        return t*t + t+t + 1.0;
    }
    // 0 otherwise
    else return 0.0;
}

float BandLimitedOsc::process() {
    float v1 = 0.f, v2 = 0.f, pointer = 0.f, numh = 0.f, pos = 0.f;
    float inc2 = 0.f, fade = 0.f, value = 0.f, maxHarms = 0.f;
    float phase = m_pointer_pos / 2*M_PI;

    switch (m_wavetype) {
        // Sine
        case 0:
            value = sinf(m_twopi * m_pointer_pos);
            break;
            // state = std::sin(phase) * amplitude;
        // Triangle
        case 1:
            maxHarms = m_srOverFour / m_freq;
            if ((m_sharp * 36.f) > maxHarms)
                numh = maxHarms / 36.f;
            else
                numh = m_sharp;
            v1 = tanf(sinf(m_twopi * m_pointer_pos)) * m_oneOverPiOverTwo;
            pointer = m_pointer_pos + 0.25f;
            if (pointer >= 1.f)
                pointer -= 1.f;
            v2 = 4.f * (0.5f - fabsf(pointer - 0.5f)) - 1.f;
            value = v1 + (v2 - v1) * numh;
            break;
        // Square
        case 2:
            maxHarms = m_srOverEight / m_freq;
            numh = m_sharp * 46.f + 4.f;
            if (numh > maxHarms)
                numh = maxHarms;
            value = atanf(numh * sinf(m_twopi * m_pointer_pos)) * m_oneOverPiOverTwo;
            break;
        // Saw
        case 3:
            maxHarms = m_srOverFour / m_freq;
            numh = m_sharp * 46.f + 4.f;
            if (numh > maxHarms)
                numh = maxHarms;
            pos = m_pointer_pos + 0.5f;
            if (pos >= 1.f)
                pos -= 1.f;
            pos = pos * 2.f - 1.f;
            value = -(pos - tanhf(numh * pos) / tanhf(numh));
            break;
        // Ramp
        case 4:
            maxHarms = m_srOverFour / m_freq;
            numh = m_sharp * 46.f + 4.f;
            if (numh > maxHarms)
                numh = maxHarms;
            pos = m_pointer_pos + 0.f;
            if (pos >= 1.f)
                pos -= 1.f;
            pos = pos * 2.f - 1.f;
            value = pos - tanhf(numh * pos) / tanhf(numh);
            break;
        // Pulse
        case 5:
            
            /*
            maxHarms = m_srOverEight / m_freq;
            numh = m_sharp * 46.f + 4.f;
            if (numh > maxHarms)
                numh = maxHarms;
            if (fmodf(numh, 2.f) == 0.f)
                numh += 1.f;
            value = tanf(powf(fabsf(sinf(m_twopi * m_pointer_pos)), numh));
            //value *= m_oneOverPiOverTwo;
            value *= m_oneOverPi;
            break;
             */
            
            
            // polyblep saw
            
        {double t = m_pointer_pos / 2*M_PI;
            value = (2.0*m_pointer_pos / 2*M_PI) - 1.0;
            value -= poly_blep(t);
            break;}
            
            
        // Bi-Pulse
        case 6:
        {
            
            /* PWM with double SAW ??*/
            double pulseWidth = 0.80f;
            //maxHarms = m_srOverEight / m_freq;
             maxHarms = m_srOverFour / m_freq;
             numh = m_sharp * 46.f + 4.f;
             if (numh > maxHarms)
                 numh = maxHarms;
             pos = m_pointer_pos + 0.5f;
             if (pos >= 1.f)
                 pos -= 1.f;
             pos = pos * 2.f - 1.f;
             value = -(pos - tanhf(numh * pos) / tanhf(numh)*2);
             double saw1 = value;

            // --- phase shift on second oscillator
            maxHarms = m_srOverFour / m_freq;
            numh = m_sharp * 46.f + 4.f;
            if (numh > maxHarms)
                numh = maxHarms;
            //pos = m_pointer_pos + 0.5f;
            pos = pos + pulseWidth;
            if (pos >= 1.f)
                pos -= 1.f;
            pos = pos * 2.f - 1.f;
            value = -(pos - tanhf(numh * pos) / tanhf(numh));
            double saw2 = value/2.0;

            double squareOut = 0.5*saw1 - 0.5*saw2;

             //--- apply DC correction
            double dcCorrection = 1.0 / pulseWidth*2;

            // --- modfiy for less than 50%
            if (pulseWidth < 0.5)
                dcCorrection = 1.0 / (1.0 - pulseWidth);

            // --- apply correction
            squareOut *= dcCorrection;

            value = squareOut;
            //value = sinf(m_twopi * m_pointer_pos);
            
    

            
            break;}
            
            /*
            
             // OLIVER BI-PULSE
            maxHarms = (int)(m_srOverEight / m_freq);
            numh = floorf(m_sharp * 46.f + 4.f);
            if (numh > maxHarms)
                numh = maxHarms;
            if (fmodf(numh, 2.f) == 0.f)
                numh += 1.f;
            value = tanf(powf(sinf(m_twopi * m_pointer_pos), numh));
            value *= m_oneOverPiOverTwo;
            break;
             */
        // SAH
        case 7:
        {
            // OLIVER BI-PULSE
//           maxHarms = (int)(m_srOverEight / m_freq);
//           numh = floorf(m_sharp * 46.f + 4.f);
//           if (numh > maxHarms)
//               numh = maxHarms;
//           if (fmodf(numh, 2.f) == 0.f)
//               numh += 1.f;
//           value = tanf(powf(sinf(m_twopi * m_pointer_pos), numh));
//           value *= m_oneOverPiOverTwo;
            
            // naive saw ? plus polyblep ? ?? MAKES A NOISE
            
            
            
//            double t = m_pointer_pos / 2*M_PI;
//            float width = 0.8;
//            value = ( m_pointer_pos > width ? 1.0 : -1.0) + (width * 2.0 - 1.0);
//            value -= poly_blep(t);
//
//
//
//            double tt = fmod(m_pointer_pos/(2*M_PI),(double)1.0);
//            if (t>0.5)
//                value =  -sin(m_pointer_pos);
//            if (t<=0.5)
//                value = (2.0*tt)-1.0;
            
            // SOUNDS GOOD !! A BAND LIMITED PULSE WAV
            pulse_width = 0.7;
            phase1 = m_pointer_pos + 0.5 * pulse_width;
            phase2 = m_pointer_pos - 0.5 * pulse_width;
            
            
            maxHarms = m_srOverFour / m_freq;
            numh = m_sharp * 46.f + 4.f;
            
            if (numh > maxHarms)
                numh = maxHarms;
            pos = phase1 + 0.5f;
            if (pos >= 1.f)
                pos -= 1.f;
            pos = pos * 2.f - 1.f;
            value1 = -(pos - tanhf(numh * pos) / tanhf(numh));
            
            pos = phase2 + 0.5f;
            if (pos >= 1.f)
                pos -= 1.f;
            pos = pos * 2.f - 1.f;
            value2 = -(pos - tanhf(numh * pos) / tanhf(numh));
            
            value = value1 - value2;
            


            //state = std::sin(phase) * amplitude;
            //double width = 0.8;

            //double state = ((phase < width) * 1.0) + ((phase >= 0.5) * (phase < 0.5 + width) * -1.0);
            //value = state;
           // double width = 0.8;
            //phase = m_twopi * m_pointer_pos;
           // double pulse = (phase > width ? 1.0 : -1.0) + (width * 2.0 - 1.0);
           // value = pulse;
            
            
            /*
            float pulseWidth = 0.8;
            
            maxHarms = m_srOverFour / m_freq;
            numh = m_sharp * 46.f + 4.f;
            if (numh > maxHarms)
                numh = maxHarms;
            pos = m_pointer_pos + 0.5f;
            if (pos >= 1.f)
                pos -= 1.f;
            pos = pos * 2.f - 1.f;
            double sawtoothSample  = -(pos - tanhf(numh * pos) / tanhf(numh));

            // --- phase shift on second oscillator
            pos = m_pointer_pos + pulseWidth;
            if (pos >= 1.f)
                pos -= 1.f;
            //pos = pos * 2.f - 1.f;

            // --- generate 2nd saw: false = do not advance the clock
            maxHarms = m_srOverFour / m_freq;
            numh = m_sharp * 46.f + 4.f;
            if (numh > maxHarms)
                numh = maxHarms;
            pos = m_pointer_pos + 0.5f;
            if (pos >= 1.f)
                pos -= 1.f;
            pos = pos * 2.f - 1.f;
            double saw2  = -(pos - tanhf(numh * pos) / tanhf(numh));

            // --- subtract = 180 out of phase
            double squareOut = 0.5*sawtoothSample - 0.5*saw2;

            // --- apply DC correction
            double dcCorrection = 1.0 / pulseWidth;

            // --- modfiy for less than 50%
            if (pulseWidth < 0.5)
                dcCorrection = 1.0 / (1.0 - pulseWidth);

            // --- apply correction
            squareOut *= dcCorrection;
            
            value = squareOut;
             
             */

            
    
            
        
            break;}

            
            /*
            Take an upramping sawtooth and its inverse, a downramping sawtooth. Adding these two waves
            with a well defined delay between 0 and period (1/f)
            results in a square wave with a duty cycle ranging from 0 to 100%.
             */
            
            /*
        {maxHarms = m_srOverFour / m_freq;
            numh = m_sharp * 46.f + 4.f;
            if (numh > maxHarms)
                numh = maxHarms;
            pos = m_pointer_pos + 0.5f;
            if (pos >= 1.f)
                pos -= 1.f;
            pos = pos * 2.f - 1.f;
            double valueOne = -(pos - tanhf(numh * pos) / tanhf(numh));
            
            pos = pos + 0.5f;
            if (pos >= 1.f)
                pos -= 1.f;
            pos = pos * 2.f - 1.f;
            double valueTwo = 1/(-(pos - tanhf(numh * pos) / tanhf(numh)));
            value = valueOne + valueTwo;
            break;}
             */
            
            
            
            
            /*
            double t = m_pointer_pos / 2*M_PI;
            double valueOne = (2.0*m_pointer_pos / 2*M_PI) - 1.0;
            valueOne -= poly_blep(t);
            
            double t2 = m_pointer_pos / 2*M_PI + 0.2;
            double valueTwo = 1/((2.0*t2) - 1.0);
            valueTwo -= poly_blep(t2);
            
            value = valueOne + valueTwo;
             */
            
            
            
            
    
            
            /*
            // polybleb pulse width ?
            double t = m_pointer_pos / 2*M_PI;
            double pulseWidth = 0.8;
            t = m_pointer_pos / 2*M_PI;
            
            double t1 = t + 0.875 + 0.25 * (pulseWidth - 0.5);
            t1 -= bitwiseOrZero(t1);

            double t2 = t + 0.375 + 0.25 * (pulseWidth - 0.5);
            t2 -= bitwiseOrZero(t2);

            // Square #1
            double y = t1 < 0.5 ? 1 : -1;

            y += poly_blep(t1) - poly_blep(t2);

            t1 += 0.5 * (1 - pulseWidth);
            t1 -= bitwiseOrZero(t1);

            t2 += 0.5 * (1 - pulseWidth);
            t2 -= bitwiseOrZero(t2);

            // Square #2
            y += t1 < 0.5 ? 1 : -1;

            y += poly_blep(t1) - poly_blep(t2);

            value =  y;
            break;
            */
            
            
            
            /*
            numh = 1.f - m_sharp;
            inc2 = 1.f / (1.f / (m_freq * m_oneOverSr) * numh);
            if (m_pointer_pos >= 1.f) {
                m_pointer_pos -= 1.f;
                m_sah_pointer_pos = 0.f;
                m_sah_last_value = m_sah_current_value;
                m_sah_current_value = (rand() / (float)RAND_MAX) * 2.f - 1.f;
            }
            if (m_sah_pointer_pos < 1.f) {
                fade = 0.5f * sinf(M_PI * (m_sah_pointer_pos + 0.5f)) + 0.5f;
                value = m_sah_current_value + (m_sah_last_value - m_sah_current_value) * fade;
                m_sah_pointer_pos += inc2;
            }
            else {
                value = m_sah_current_value;
            }
            m_pointer_pos += m_freq * m_oneOverSr;
            break;
             */
             
        default:
            value = 0.f;
            break;
             
    }

    if (m_wavetype < 8) {
        m_pointer_pos += m_freq * m_oneOverSr;
        m_pointer_pos = _clip(m_pointer_pos);
    }

    return value;
}

