<?xml version="1.0" encoding="UTF-8"?>
<ui version="4.0">
 <class>Dialog</class>
 <widget class="QDialog" name="Dialog">
  <property name="geometry">
   <rect>
    <x>0</x>
    <y>0</y>
    <width>400</width>
    <height>300</height>
   </rect>
  </property>
  <property name="windowTitle">
   <string>Dialog</string>
  </property>
  <widget class="QCheckBox" name="cNumeric">
   <property name="geometry">
    <rect>
     <x>300</x>
     <y>40</y>
     <width>89</width>
     <height>20</height>
    </rect>
   </property>
   <property name="text">
    <string>Numeric</string>
   </property>
  </widget>
  <widget class="QPushButton" name="bCancel">
   <property name="geometry">
    <rect>
     <x>270</x>
     <y>220</y>
     <width>121</width>
     <height>61</height>
    </rect>
   </property>
   <property name="text">
    <string>Cancel</string>
   </property>
  </widget>
  <widget class="QPushButton" name="bFieldValues">
   <property name="geometry">
    <rect>
     <x>160</x>
     <y>60</y>
     <width>71</width>
     <height>32</height>
    </rect>
   </property>
   <property name="text">
    <string>Values</string>
   </property>
  </widget>
  <widget class="QPushButton" name="bOK">
   <property name="geometry">
    <rect>
     <x>10</x>
     <y>220</y>
     <width>121</width>
     <height>61</height>
    </rect>
   </property>
   <property name="text">
    <string>Plot</string>
   </property>
  </widget>
  <widget class="QLabel" name="label">
   <property name="geometry">
    <rect>
     <x>60</x>
     <y>40</y>
     <width>41</width>
     <height>21</height>
    </rect>
   </property>
   <property name="text">
    <string>Field:</string>
   </property>
  </widget>
  <widget class="QLabel" name="lExpName">
   <property name="geometry">
    <rect>
     <x>180</x>
     <y>0</y>
     <width>59</width>
     <height>16</height>
    </rect>
   </property>
   <property name="text">
    <string>ExpName</string>
   </property>
  </widget>
  <widget class="QComboBox" name="cField">
   <property name="geometry">
    <rect>
     <x>100</x>
     <y>40</y>
     <width>191</width>
     <height>26</height>
    </rect>
   </property>
  </widget>
  <widget class="QCheckBox" name="cLines">
   <property name="geometry">
    <rect>
     <x>30</x>
     <y>110</y>
     <width>89</width>
     <height>20</height>
    </rect>
   </property>
   <property name="text">
    <string>Draw lines</string>
   </property>
   <property name="checked">
    <bool>true</bool>
   </property>
   <property name="tristate">
    <bool>false</bool>
   </property>
  </widget>
  <widget class="QLabel" name="label_2">
   <property name="geometry">
    <rect>
     <x>30</x>
     <y>160</y>
     <width>71</width>
     <height>20</height>
    </rect>
   </property>
   <property name="text">
    <string>Min. reads:</string>
   </property>
  </widget>
  <widget class="QSpinBox" name="sMinReads">
   <property name="geometry">
    <rect>
     <x>100</x>
     <y>160</y>
     <width>47</width>
     <height>24</height>
    </rect>
   </property>
   <property name="maximum">
    <number>1000</number>
   </property>
   <property name="value">
    <number>4</number>
   </property>
  </widget>
  <widget class="QCheckBox" name="cSort">
   <property name="geometry">
    <rect>
     <x>60</x>
     <y>20</y>
     <width>89</width>
     <height>20</height>
    </rect>
   </property>
   <property name="text">
    <string>Sort</string>
   </property>
   <property name="checked">
    <bool>true</bool>
   </property>
  </widget>
  <widget class="QPushButton" name="bMetaData">
   <property name="geometry">
    <rect>
     <x>270</x>
     <y>130</y>
     <width>115</width>
     <height>32</height>
    </rect>
   </property>
   <property name="text">
    <string>MetaData</string>
   </property>
  </widget>
 </widget>
 <resources/>
 <connections/>
</ui>
