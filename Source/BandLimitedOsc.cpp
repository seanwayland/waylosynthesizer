

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
}

void BandLimitedOsc::setup(float sampleRate) {
    m_sampleRate = sampleRate;
    m_oneOverSr = 1.f / m_sampleRate;
    m_twopi = 2.f * M_PI;
    m_oneOverPiOverTwo = 1.f / (M_PI / 2.f);
    m_oneOverPi = 1.f / M_PI;
    m_srOverFour = m_sampleRate / 4.f;
    m_srOverEight = m_sampleRate / 8.f;
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

    switch (m_wavetype) {
        // Sine
        case 0:
            value = sinf(m_twopi * m_pointer_pos);
            break;
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
             
        // Bi-Pulse
        case 6:
        {
            
            // polyblep saw
            double t = m_pointer_pos / 2*M_PI;
            value = (2.0*m_pointer_pos / 2*M_PI) - 1.0;
            value -= poly_blep(t);
            break;}
            /*
            
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
        default:
            value = 0.f;
            break;
             */
    }

    if (m_wavetype < 7) {
        m_pointer_pos += m_freq * m_oneOverSr;
        m_pointer_pos = _clip(m_pointer_pos);
    }

    return value;
}

