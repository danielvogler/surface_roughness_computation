# find files for b.stl side and rename to a.stl

num=0
for i in `find -name "*b.stl"`; do
    let num=$((num+1))
    cp $i $i.copy
done

for i in `find -name "*copy"`; do
    rename 's/b.stl.copy/a.stl/' $i
done
