from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import re
import time
from datetime import datetime
import os
import math
from sklearn.metrics import roc_auc_score
from sklearn.metrics import confusion_matrix

import tensorflow as tf
import numpy as np
import tensorflow.contrib.slim as slim


from utils import argparser, logging
from dataset import DataSet
from model import get_placeholders, inference



def test( FLAGS ):
  # read data
  dataset = DataSet( fpath = FLAGS.test_file, 
                      seqlen = FLAGS.seq_len,
                      n_classes = FLAGS.num_classes,
                      num_feature = FLAGS.num_feature,
                      is_raw = FLAGS.is_raw,
                      need_shuffle = False )
  # set character set size
  FLAGS.charset_size = dataset.charset_size

  with tf.Graph().as_default():
    # placeholder
    placeholders = get_placeholders(FLAGS)
    
    # get inference
    pred, layers = inference( placeholders['data'], FLAGS, 
                      for_training=False )

    prob = tf.nn.softmax(pred)
    # calculate prediction
    _hit_op = tf.equal( tf.argmax(pred, 1), tf.argmax(placeholders['labels'], 1))
    hit_op = tf.reduce_sum( tf.cast( _hit_op ,tf.float32 ) )

    # create saver
    saver = tf.train.Saver()

    # summary 
    summary_op = tf.summary.merge_all()

    with tf.Session() as sess:
      # load model
      ckpt = tf.train.latest_checkpoint( os.path.dirname( FLAGS.checkpoint_path ) )
      if tf.train.checkpoint_exists( ckpt ):
        saver.restore( sess, ckpt )
        global_step = ckpt.split('/')[-1].split('-')[-1]
        logging('Succesfully loaded model from %s at step=%s.' %
              (ckpt, global_step), FLAGS)
      else:
        logging("[ERROR] Checkpoint not exist", FLAGS)
        return

      # iter batch
      hit_count = 0.0
      total_count = 0
      pred_list = []
      label_list = []

      logging("%s: starting test." % (datetime.now()), FLAGS)
      start_time = time.time()
      total_batch_size = math.ceil( dataset._num_data / FLAGS.batch_size )

      fout = open(FLAGS.out_file,'w')

      for step, (data, labels) in enumerate(dataset.iter_once( FLAGS.batch_size )):
        hits, pred_val = sess.run( [hit_op, prob], feed_dict={
          placeholders['data']: data,
          placeholders['labels']: labels
        })

        hit_count += np.sum( hits )
        total_count += len( data )

        for i, p in enumerate(pred_val):
          pred_list.append( p[0] )
          label_list.append( labels[i][0] )

        if step % FLAGS.log_interval == 0:
          duration = time.time() - start_time
          sec_per_batch = duration / FLAGS.log_interval
          examples_per_sec = FLAGS.batch_size / sec_per_batch
          logging('%s: [%d batches out of %d] (%.1f examples/sec; %.3f'
                'sec/batch)' % (datetime.now(), step, total_batch_size,
                                examples_per_sec, sec_per_batch), FLAGS)
          start_time = time.time()


      # micro precision
      # logging("%s: micro-precision = %.5f" % 
      #       (datetime.now(), (hit_count/total_count)), FLAGS)

      auc_val = roc_auc_score(label_list, pred_list)
      logging("%s: micro-precision = %.5f, auc = %.5f" % 
            (datetime.now(), (hit_count/total_count), auc_val), FLAGS)

      pred_y = [1 if i > FLAGS.threshold else 0 for i in pred_list]

      TN, FP, FN, TP = confusion_matrix(label_list, pred_y, labels=[0, 1]).ravel()

      Sensitivity = round((TP/(TP+FN)), 4) if TP+FN > 0 else 0
      Specificity = round(TN/(FP+TN), 4) if FP+TN > 0 else 0
      Precision = round(TP/(TP+FP), 4) if TP+FP > 0 else 0
      Accuracy = round((TP+TN)/(TP+FP+TN+FN), 4)
      MCC = round(((TP*TN)-(FP*FN))/(math.sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))),
                  4) if TP+FP > 0 and FP+TN > 0 and TP+FN and TN+FN else 0
      F1 = round((2*TP)/((2*TP)+FP+FN), 4)

      fout.write(f"TP={TP}, FP={FP}, TN={TN}, FN={FN}, Sens={Sensitivity}, Spec={Specificity}, Prec={Precision}, Acc={Accuracy}, MCC={MCC}, F1={F1}\n")

      logging(f"TP={TP}, FP={FP}, TN={TN}, FN={FN}, Sens={Sensitivity}, Spec={Specificity}, Prec={Precision}, Acc={Accuracy}, MCC={MCC}, F1={F1}", FLAGS)




if __name__ == '__main__':
  FLAGS = argparser()
  FLAGS.is_training = False
  logging(str(FLAGS), FLAGS)
  test( FLAGS )
