#include <SoftwareSerial.h>
SoftwareSerial BTSerial(10, 11);  // RX, TX

const int ecgPin = A0;

void setup() {
  BTSerial.begin(9600);  // HC-05 default baud
  Serial.begin(9600);    // Optional: for USB Serial Monitor
}

void loop() {
  int ecgValue = analogRead(ecgPin);
  BTSerial.println(ecgValue);	// Send to HC-05
  Serial.println(ecgValue);	// Optional: send to Serial Monitor
  delay(3);
}

