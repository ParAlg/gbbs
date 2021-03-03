#pragma once

#include "gbbs/gbbs.h"

namespace gbbs {

template<class E>
struct BBQueueNode {
  E data;
  struct BBQueueNode *next = nullptr;
  struct BBQueueNode *prev = nullptr;
  BBQueueNode(){}
  BBQueueNode(E _data):data(_data){}
  void erase(){
    if(next && prev){
        next->prev, prev->next = prev, next;
    }else if(next){
        next->prev = nullptr;
    }else{
        prev->next = nullptr;
    }
      
  }
};

template<class E>
struct BBQueue{
    BBQueueNode<E> *head = nullptr;
    BBQueueNode<E> *tail = nullptr;
    size_t n = 0;
    size_t size(){return n;}
    bool empty(){return n==0;}


    BBQueue(){}

    E front(){
        return head->data;
    }

    E pop_front(){
        // if(n==0) return; can't be empty
        BBQueueNode<E> *tmp = head->next;
        delete head;
        head = tmp;
        n--;
    }

    // assume node is in the queue
    void erase(BBQueueNode<E> *node){
        node->erase();
        delete node;
        n--;
    }

    bool erase(E data){
        BBQueueNode<E>* node = head;
        while(node){
            if(node->data == data){
                erase(node);
                return true;
            }
        }
        return false;
    }

    void push_back(E data){
        BBQueueNode<E>* node = new BBQueueNode<E>(data);
        if(n==0){
            head = node;
            tail = node;
        }else{
            node->prev = tail;
            tail->next = node;
            tail = node;
        }
        n++;
    }
};

}//end of gbbs