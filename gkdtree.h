// MIT License
// Copyright (c) 2018 Bo Yang
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to
// deal in the Software without restriction, including without limitation the
// rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
// sell copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.

#ifndef __GENERIC_KD_TREE_DEFINITION_H__
#define __GENERIC_KD_TREE_DEFINITION_H__
#include <algorithm>
#include <iterator>
#include <memory>
#include <tuple>
#include <type_traits>

namespace GenericKDTree {

//
// default values for all types
//
namespace type_traits {

// use int64_t as default type for coordinate
template <typename Tp> struct CoordinateType { typedef int64_t type; };

// default : 2D problem
template <typename Tp> struct Dimension { constexpr static size_t value = 2; };

// default: max tree depth
template <typename Tp> struct MaxLevel { constexpr static size_t value = 64; };

// default: max num of elements associated to a tree node
template <typename Tp> struct NodeSize { constexpr static size_t value = 10; };
}

// must be specialized for T & Axis, otherwise compile error:
template <typename T, size_t Axis>
typename type_traits::CoordinateType<T>::type GetLowerCoordinate(const T &p);

template <typename T, size_t Axis>
typename type_traits::CoordinateType<T>::type GetUpperCoordinate(const T &p);

namespace DimensionalOpImpl {

// helper function to check if current axis
// is along the last dimension
template <typename T, size_t CurrAxis> constexpr bool IsLastDimension() {

  return (CurrAxis == type_traits::Dimension<T>::value - 1);
}

// helper types to iterate over all dimensions
template <typename T, size_t CurrAxis, bool> struct NextAxisHelper {

  constexpr static size_t value = CurrAxis + 1;
};

// reached max dimension, reset to the first dimension
template <typename T, size_t CurrAxis>
struct NextAxisHelper<T, CurrAxis, true> {

  constexpr static size_t value = 0;
};

template <typename T, size_t CurrAxis> constexpr size_t NextAxis() {

  static_assert(CurrAxis < type_traits::Dimension<T>::value,
                "Wrong object dimension!");

  return NextAxisHelper<T, CurrAxis, IsLastDimension<T, CurrAxis>()>::value;
}
}

namespace DefaultOverlapCheckImpl {

// default helper types to check overlap on all dimensions

// simple version, assuming rectangular objects.
// should be fucntion that check overlap along a given
// axis (dimension)
template <typename T, size_t CurrAxis>
bool DefaultDimensionalOverlapCheckHelper(const T &a, const T &b) {
  return !(
      GetLowerCoordinate<T, CurrAxis>(a) > GetUpperCoordinate<T, CurrAxis>(b) ||
      GetLowerCoordinate<T, CurrAxis>(b) > GetUpperCoordinate<T, CurrAxis>(a));
}

//
// Recursive template to check overlap on all dimensions
//
template <typename T, size_t CurrAxis, bool LastDimen>
struct DefaultOverlapChecker {

  static bool apply(const T &a, const T &b) {

    typedef DefaultOverlapChecker<
        T, DimensionalOpImpl::NextAxis<T, CurrAxis>(),
        DimensionalOpImpl::IsLastDimension<
            T, DimensionalOpImpl::NextAxis<T, CurrAxis>()>()>
        DefaultNextDimensionOverlapChecker;

    // curr dimension overlapped and next dimension overlapped
    return OverlappedOnCurrentDimension(a, b) &&
           DefaultNextDimensionOverlapChecker::apply(a, b);
  }

private:
  static bool OverlappedOnCurrentDimension(const T &a, const T &b) {
    return DefaultDimensionalOverlapCheckHelper<T, CurrAxis>(a, b);
  }
};

// check along the last dimesion:
template <typename T, size_t CurrAxis>
struct DefaultOverlapChecker<T, CurrAxis, true> {

  static bool apply(const T &a, const T &b) {
    return DefaultDimensionalOverlapCheckHelper<T, CurrAxis>(a, b);
  }
};
}

// Check overlap with default impl.
template <typename T> bool Overlapped(const T &a, const T &b) {

  // start from the first dimension (0)
  // template will be recursively expanded to apply checks on
  // all dimensions.
  return DefaultOverlapCheckImpl::DefaultOverlapChecker<
      T, 0, DimensionalOpImpl::IsLastDimension<T, 0>()>::apply(a, b);
}

template <typename T, size_t Axis, typename InputIterator>
std::tuple<typename type_traits::CoordinateType<T>::type,
           typename type_traits::CoordinateType<T>::type>
CoordinateRange(InputIterator begin, InputIterator end) {

  typedef typename type_traits::CoordinateType<T>::type CoordinateType;

  CoordinateType min = GetLowerCoordinate<T, Axis>(
      *(std::min_element(begin, end, [](const T &a, const T &b) {
        return GetLowerCoordinate<T, Axis>(a) < GetLowerCoordinate<T, Axis>(b);
      })));

  CoordinateType max = GetUpperCoordinate<T, Axis>(
      *(std::max_element(begin, end, [](const T &a, const T &b) {
        return GetUpperCoordinate<T, Axis>(a) < GetUpperCoordinate<T, Axis>(b);
      })));

  return std::make_tuple(min, max);
}

template <typename NodeType, size_t CurrAxis, typename IteratorType>
struct TreeType {

  typedef TreeType<NodeType, DimensionalOpImpl::NextAxis<NodeType, CurrAxis>(),
                   IteratorType>
      ChildNodeType;

  typedef typename type_traits::CoordinateType<NodeType>::type CoordinateType;
  typedef IteratorType RangeIteratorType;

  TreeType() = delete;

  TreeType(RangeIteratorType begin, RangeIteratorType end, size_t level = 0) {
    BuildIndex(begin, end, level);
  }

  TreeType(TreeType &&t)
      : lower_tree_(std::move(t.lower_tree_)),
        upper_tree_(std::move(t.upper_tree_)), data_(std::move(data_)),
        span_(std::move(span_)) {}

  template <typename ContainerT>
  ContainerT QueryOverlappedObjects(const NodeType &window) const {

    if (GetLowerCoordinate<NodeType, CurrAxis>(window) > std::get<1>(span_)) {
      return ContainerT();
    }

    if (GetUpperCoordinate<NodeType, CurrAxis>(window) < std::get<0>(span_)) {
      return ContainerT();
    }

    ContainerT results;

    if (lower_tree_) {
      results = lower_tree_->QueryOverlappedObjects(window);
    }

    if (upper_tree_) {
      ContainerT tmp = upper_tree_->QueryOverlappedObjects(window);
      results.insert(results.end(), tmp.begin(), tmp.end());
    }

    ContainerT tmp = OverlappedWithCurrentDataSet<ContainerT>(window);
    results.insert(results.end(), tmp.begin(), tmp.end());

    return results;
  }

  template <typename Predicate>
  void ForEachObjectOverlap(const NodeType &window, Predicate Pred) {

    if (lower_tree_)
      lower_tree_->ForEachObjectOverlap(window, Pred);
    ForCurrentDataSetOverlappedWith(window, Pred);
    if (upper_tree_)
      upper_tree_->ForEachObjectOverlap(window, Pred);
  }

private:
  template <typename ContainerT>
  ContainerT OverlappedWithCurrentDataSet(const NodeType &window) const {

    ContainerT results;
    auto pred = [&results](const NodeType &obj) { results.push_back(obj); };

    ForCurrentDataSetOverlappedWith(window, pred);

    return results;
  }

  template <typename Predicate>
  void ForCurrentDataSetOverlappedWith(const NodeType &window, Predicate Pred) {

    // quick check if there is overlap:
    std::tuple<CoordinateType, CoordinateType> data_span_ =
        CoordinateRange<NodeType, CurrAxis>(std::get<0>(data_),
                                            std::get<1>(data_));

    if (GetLowerCoordinate<NodeType, CurrAxis>(window) >
            std::get<1>(data_span_) ||
        GetUpperCoordinate<NodeType, CurrAxis>(window) <
            std::get<0>(data_span_)) {

      return;
    }

    std::for_each(std::get<0>(data_), std::get<1>(data_),
                  [&window, &Pred](NodeType &obj) {
                    if (Overlapped(window, obj))
                      Pred(obj);
                  });
  }

  bool BuildIndex(RangeIteratorType begin, RangeIteratorType end,
                  size_t level) {

    /* update coordinates span in all dimensions */

    span_ = CoordinateRange<NodeType, CurrAxis>(begin, end);

    /* reached max level, stop partitoning. */
    if (level == type_traits::MaxLevel<NodeType>::value ||
        std::distance(begin, end) < type_traits::NodeSize<NodeType>::value) {

      std::get<0>(data_) = begin;
      std::get<1>(data_) = end;

      return true;
    }

    CoordinateType partition_line =
        (std::get<0>(span_) + std::get<1>(span_)) / 2;

    //
    // absolutely below the partiton_line
    // [begin .... | below_pos ... end )
    //             ^ partition_line
    //
    auto below_pos =
        std::partition(begin, end, [partition_line](const NodeType &n) {
          return GetUpperCoordinate<NodeType, CurrAxis>(n) < partition_line;
        });

    // build lower sub-tree
    if (below_pos != end) {

      lower_tree_ = std::unique_ptr<ChildNodeType>(
          new ChildNodeType(begin, below_pos, level + 1));
    }

    // absolutely above the partition_line:
    auto upper_pos =
        std::partition(below_pos, end, [partition_line](const NodeType &n) {
          return GetLowerCoordinate<NodeType, CurrAxis>(n) <= partition_line;
        });

    if (upper_pos != end) {
      upper_tree_ = std::unique_ptr<ChildNodeType>(
          new ChildNodeType(upper_pos, end, level + 1));
    }

    // [below_pos, upper_pos]
    std::get<0>(data_) = below_pos;
    std::get<1>(data_) = upper_pos;
    return true;
  }

private:
  // data members
  std::unique_ptr<ChildNodeType> lower_tree_;
  std::unique_ptr<ChildNodeType> upper_tree_;
  std::tuple<RangeIteratorType, RangeIteratorType> data_;
  std::tuple<CoordinateType, CoordinateType> span_;
};

// helper function for making a kd-tree:
template <typename IteratorType>
TreeType<typename std::iterator_traits<IteratorType>::value_type, 0,
         IteratorType>
make_tree(IteratorType begin, IteratorType end) {

  static_assert(
      std::is_same<
          typename std::iterator_traits<IteratorType>::iterator_category,
          std::random_access_iterator_tag>::value,
      "Must be random access iterator!");

  typedef typename std::iterator_traits<IteratorType>::value_type NodeType;
  return TreeType<NodeType, 0, IteratorType>(begin, end);
}
}

#endif
