SUBROUTINE pack_left_right_buffers

  ! These array modifications still need to be added on, plus the donor data location changes as in update_halo
  IF(field_type.EQ.CELL_DATA) THEN
    x_inc=0
    y_inc=0
  ENDIF
  IF(field_type.EQ.VERTEX_DATA) THEN
    x_inc=1
    y_inc=1
  ENDIF
  IF(field_type.EQ.X_FACE_DATA) THEN
    x_inc=1
    y_inc=0
  ENDIF
  IF(field_type.EQ.Y_FACE_DATA) THEN
    x_inc=0
    y_inc=1
  ENDIF

  ! Pack real data into buffers
  IF(parallel%task.EQ.chunks(chunk)%task) THEN
    size=(1+(chunks(chunk)%field%y_max+y_inc+depth)-(chunks(chunk)%field%y_min-depth))*depth
!$OMP PARALLEL
    IF(chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) THEN
!$OMP DO PRIVATE(index)
      DO k=chunks(chunk)%field%y_min-depth,chunks(chunk)%field%y_max+y_inc+depth
        DO j=1,depth
          index=j+(k+depth-1)*depth
          left_snd_buffer(index)=field(chunks(chunk)%field%x_min+x_inc-1+j,k)
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
    IF(chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) THEN
!$OMP DO PRIVATE(index)
      DO k=chunks(chunk)%field%y_min-depth,chunks(chunk)%field%y_max+y_inc+depth
        DO j=1,depth
          index=j+(k+depth-1)*depth
          right_snd_buffer(index)=field(chunks(chunk)%field%x_max+1-j,k)
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
!$OMP END PARALLEL

END SUBROUTINE pack_left_right_buffers

SUBROUTINE unpack_left_right_buffers

 ! Unpack buffers in halo cells
  IF(parallel%task.EQ.chunks(chunk)%task) THEN
!$OMP PARALLEL
    IF(chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) THEN
!$OMP DO PRIVATE(index)
      DO k=chunks(chunk)%field%y_min-depth,chunks(chunk)%field%y_max+y_inc+depth
        DO j=1,depth
          index=j+(k+depth-1)*depth
          field(chunks(chunk)%field%x_min-j,k)=left_rcv_buffer(index)
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
    IF(chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) THEN
!$OMP DO PRIVATE(index)
      DO k=chunks(chunk)%field%y_min-depth,chunks(chunk)%field%y_max+y_inc+depth
        DO j=1,depth
          index=j+(k+depth-1)*depth
          field(chunks(chunk)%field%x_max+x_inc+j,k)=right_rcv_buffer(index)
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
!$OMP END PARALLEL
  ENDIF

END SUBROUTINE unpack_left_right_buffers

SUBROUTINE pack_top_bottom_buffers

  IF(parallel%task.EQ.chunks(chunk)%task) THEN
    size=(1+(chunks(chunk)%field%x_max+x_inc+depth)-(chunks(chunk)%field%x_min-depth))*depth
!$OMP PARALLEL
    IF(chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) THEN
      DO k=1,depth
!$OMP DO PRIVATE(index)
        DO j=chunks(chunk)%field%x_min-depth,chunks(chunk)%field%x_max+x_inc+depth
          index=j+depth+(k-1)*(chunks(chunk)%field%x_max+x_inc+(2*depth))
          bottom_snd_buffer(index)=field(j,chunks(chunk)%field%y_min+y_inc-1+k)
        ENDDO
!$OMP END DO
      ENDDO
    ENDIF
    IF(chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) THEN
      DO k=1,depth
!$OMP DO PRIVATE(index)
        DO j=chunks(chunk)%field%x_min-depth,chunks(chunk)%field%x_max+x_inc+depth
          index=j+depth+(k-1)*(chunks(chunk)%field%x_max+x_inc+(2*depth))
          top_snd_buffer(index)=field(j,chunks(chunk)%field%y_max+1-k)
        ENDDO
!$OMP END DO
      ENDDO
    ENDIF
!$OMP END PARALLEL

END SUBROUTINE pack_top_bottom_buffers

SUBROUTINE unpack_top_bottom_buffers

  ! Unpack buffers in halo cells
  IF(parallel%task.EQ.chunks(chunk)%task) THEN
!$OMP PARALLEL
    IF(chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) THEN
      DO k=1,depth
!$OMP DO PRIVATE(index)
        DO j=chunks(chunk)%field%x_min-depth,chunks(chunk)%field%x_max+x_inc+depth
          index=j+depth+(k-1)*(chunks(chunk)%field%x_max+x_inc+(2*depth))
          field(j,chunks(chunk)%field%y_min-k)=bottom_rcv_buffer(index)
        ENDDO
!$OMP END DO
      ENDDO
    ENDIF
    IF(chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) THEN
      DO k=1,depth
!$OMP DO PRIVATE(index)
        DO j=chunks(chunk)%field%x_min-depth,chunks(chunk)%field%x_max+x_inc+depth
          index=j+depth+(k-1)*(chunks(chunk)%field%x_max+x_inc+(2*depth))
          field(j,chunks(chunk)%field%y_max+y_inc+k)=top_rcv_buffer(index)
        ENDDO
!$OMP END DO
      ENDDO
    ENDIF
!$OMP END PARALLEL
  ENDIF

END SUBROUTINE unpack_top_bottom_buffers
