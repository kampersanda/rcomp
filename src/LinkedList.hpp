/**
 * @file LinkedList.hpp
 * @author Shunsuke Kanda (kampersanda)
 */
#pragma once

#include "abort_if.hpp"
#include "basics.hpp"

namespace rcomp {

/**
 * A class for bi-directional linked list.
 *
 * @tparam T The data type.
 */
template <class T>
class LinkedList {
  public:
    using data_type = T;

    class node_type {
      public:
        using data_type = T;

      private:
        data_type  m_data;
        node_type* m_next = nullptr;
        node_type* m_prev = nullptr;

      public:
        node_type() = default;

        size_type get_memory_in_bytes() const {
            return m_data.get_memory_in_bytes() + sizeof(m_next) + sizeof(m_prev);
        }

        void set_data(data_type&& v) {
            m_data = std::move(v);
        }
        void set_next(node_type* v) {
            m_next = v;
        }
        void set_prev(node_type* v) {
            m_prev = v;
        }

        data_type& get_data() {
            return m_data;
        }
        const data_type& get_data() const {
            return m_data;
        }
        node_type* get_next() const {
            return m_next;
        }
        node_type* get_prev() const {
            return m_prev;
        }
    };

  private:
    node_type* m_head = nullptr;

  public:
    LinkedList() = default;

    virtual ~LinkedList() {
        if (m_head == nullptr) {
            return;
        }

        node_type* x = m_head->get_next();
        while (x != m_head) {
            x = remove(x);
        }

        delete m_head;
        m_head = nullptr;
    }

    //! Copy constructor (deleted)
    LinkedList(const LinkedList&) = delete;

    //! Copy constructor (deleted)
    LinkedList& operator=(const LinkedList&) = delete;

    //! Move constructor
    LinkedList(LinkedList&&) noexcept = default;

    //! Move constructor
    LinkedList& operator=(LinkedList&&) noexcept = default;

    void clear() {
        m_head = new node_type();
        m_head->set_next(m_head);
        m_head->set_prev(m_head);
    }

    node_type* get_head() const {
        return m_head;
    }

    // FROM x1 -> x2 TO x1 -> x3 -> x2
    // RETURN x3
    node_type* insert_after(node_type* x1) {
        node_type* x3 = new node_type();
        node_type* x2 = x1->get_next();
        x3->set_prev(x1);
        x3->set_next(x2);
        x1->set_next(x3);
        x2->set_prev(x3);
        return x3;
    }

    // FROM x1 -> x2 -> x3 TO x1 -> x3
    // RETURN x3
    node_type* remove(node_type* x2) {
        ABORT_IF(x2 == m_head);
        node_type* x1 = x2->get_prev();
        node_type* x3 = x2->get_next();
        x1->set_next(x3);
        x3->set_prev(x1);
        delete x2;
        return x3;
    }

    // not in O(1) time
    size_type get_memory_in_bytes() const {
        if (m_head == nullptr) {
            return 0;
        }
        size_type mem = m_head->get_memory_in_bytes();
        for (node_type* x = m_head->get_next(); x != m_head; x = x->get_next()) {
            mem += x->get_memory_in_bytes();
        }
        return mem;
    }
};

}  // namespace rcomp
