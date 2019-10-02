/**
 * @file define.h
 * @brief Define the flags
 * @author Keren Zhu
 * @date 09/30/2019
 */

#ifndef IDEAPLACE_DEFINE_H_
#define IDEAPLACE_DEFINE_H_

//#define NODEBUG

#ifdef NODEBUG
#define AT(vec, idx) vec[idx]
#else
#define AT(vec, idx) vec.at(idx)
#endif /// NODEBUG

#endif /// IDEAPLACE_DEFINE_H
