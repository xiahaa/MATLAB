function qpqc = prod_qpqconj(q1,q2)
    qp  = q2m_left(q1)*q2;
    qpqc = q2m_left(qp)*conjugateq(q1);
end