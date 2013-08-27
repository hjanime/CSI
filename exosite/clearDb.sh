rm -r exo/migrations/
psql exosite -f test.sql
./manage.py schemamigration exo --initial
./manage.py migrate exo --delete-ghost-migrations
./manage.py schemamigration exo --initial
./manage.py migrate exo

