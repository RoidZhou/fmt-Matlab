## fmt-Matlab
一、原理
将笛卡尔坐标系下的旋转和缩放转化为新坐标系下的平移，通过相位相关求得平移量就得到了缩放倍率和旋转角度。根据倍率和旋转角度做矫正，再直接相位相关求得平移量。于是就得到了两幅图像的相对位移、旋转和缩放，可以用于图像配准。

二、步骤
1. 产生两个等大的正方形图像块

2. 对这两个图像块分别做傅里叶变换、对数极坐标变换（统称傅里叶梅林变换），下图就是变换结果

![](https://github.com/RoidZhou/fmt-Matlab/blob/main/img/fmtlp.jpg)

3. 再做傅里叶变换，然后相位相关，再逆变换就得到了响应图

4. 寻找响应图最大值位置，然后查表得到旋转角度和缩放倍率

   以上部分求得旋转、缩放，还有平移没有求

   以下部分没有写代码

5. 同样适用相位相关求平移
   下图是计算出的响应图像

![](https://github.com/RoidZhou/fmt-Matlab/blob/main/img/fmtrsp.jpg)
